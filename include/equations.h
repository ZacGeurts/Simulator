#ifndef EQUATIONS_H
#define EQUATIONS_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <functional>
#include <iomanip>
#include <memory>
#include <omp.h>
#include <fstream>
#include <mutex>
#include <atomic>
#include <chrono>
#include <random>
#include <algorithm>
#include <limits>
#include <ctime>
#include <iostream>

// Thread-safe logger
class Logger {
private:
    std::mutex mutex;
    std::vector<std::string> buffer;
    const size_t max_buffer_size = 100;

public:
    void log(const std::string& message) {
        auto now = std::time(nullptr);
        std::stringstream ss;
        ss << "[" << std::put_time(std::localtime(&now), "%Y-%m-%d %H:%M:%S") << "] " << message;

        std::lock_guard<std::mutex> lock(mutex);
        buffer.push_back(ss.str());

        if (buffer.size() >= max_buffer_size) {
            flush();
        }
    }

    void flush() {
        std::lock_guard<std::mutex> lock(mutex);
        if (buffer.empty()) return;

        std::ofstream log_file("simulation.log", std::ios::app);
        if (!log_file.is_open()) {
            std::cerr << "Failed to open simulation.log" << std::endl;
            return;
        }
        for (const auto& msg : buffer) {
            log_file << msg << std::endl;
            std::cout << msg << std::endl;
        }
        buffer.clear();
    }
};

// Tensor implementation
class Tensor {
public:
    Tensor() {
        std::fill(&data[0][0], &data[0][0] + 9, 0.0L);
    }

    void set_diagonal(long double value) {
        data[0][0] = data[1][1] = data[2][2] = value;
    }

    long double data[3][3];
};

// Particle implementation
class Particle {
public:
    Particle() : x(0.0L), y(0.0L), z(0.0L), vx(0.0L), vy(0.0L), vz(0.0L), ax(0.0L), ay(0.0L), az(0.0L) {}

    long double x, y, z;
    long double vx, vy, vz;
    long double ax, ay, az;
};

// Equation base class
class Equation {
public:
    virtual ~Equation() = default;
    virtual void compute(class Simulation& sim) = 0;
    virtual std::string output(const class Simulation& sim) const = 0;
    virtual std::string name() const = 0;
    virtual bool is_valid() const = 0;
    virtual void setVisScale(long double s) = 0;
    virtual long double getVisScale() const = 0;
    virtual void setVisColorIntensity(long double c) = 0;
    virtual long double getVisColorIntensity() const = 0;
};

// Forward declaration of Simulation
class Simulation;

// ExpressionParser class
class ExpressionParser {
private:
    std::map<std::string, std::function<long double(long double)>> functions;
    std::map<std::string, std::function<long double(long double, long double)>> binary_functions;
    std::map<std::string, long double> vars;
    std::string expression;
    size_t pos;

    std::string preprocess(const std::string& expr) {
        std::string result;
        bool last_was_digit = false;
        bool last_was_var = false;
        for (char c : expr) {
            if (c == '^') {
                result += "**";
                last_was_digit = last_was_var = false;
                continue;
            }
            if (std::isspace(c)) continue;
            if (std::isalpha(c) && last_was_digit) {
                result += "*";
            } else if (std::isdigit(c) && last_was_var) {
                result += "*";
            }
            result += c;
            last_was_digit = std::isdigit(c) || c == '.' || c == 'e';
            last_was_var = std::isalpha(c);
        }
        return result;
    }

    void skipWhitespace() {
        while (pos < expression.length() && std::isspace(expression[pos])) {
            ++pos;
        }
    }

    long double parseExpression() {
        long double result = parseTerm();
        while (pos < expression.length()) {
            skipWhitespace();
            char op = expression[pos];
            if (op != '+' && op != '-') break;
            ++pos;
            long double term = parseTerm();
            result = (op == '+') ? result + term : result - term;
            if (std::isinf(result) || std::isnan(result)) {
                throw std::runtime_error("Numerical overflow in expression");
            }
        }
        return result;
    }

    long double parseTerm() {
        long double result = parseFactor();
        while (pos < expression.length()) {
            skipWhitespace();
            char op = expression[pos];
            if (op == '*' && pos + 1 < expression.length() && expression[pos + 1] == '*') {
                pos += 2;
                long double factor = parseFactor();
                result = std::pow(result, factor);
            } else if (op != '*' && op != '/') {
                break;
            } else {
                ++pos;
                long double factor = parseFactor();
                if (op == '/' && factor == 0.0L) {
                    throw std::runtime_error("Division by zero");
                }
                result = (op == '*') ? result * factor : result / factor;
            }
            if (std::isinf(result) || std::isnan(result)) {
                throw std::runtime_error("Numerical overflow in term");
            }
        }
        return result;
    }

    long double parseFactor() {
        skipWhitespace();
        if (pos >= expression.length()) {
            throw std::runtime_error("Unexpected end of expression");
        }

        char c = expression[pos];
        if (c == '(') {
            ++pos;
            long double result = parseExpression();
            skipWhitespace();
            if (pos >= expression.length() || expression[pos] != ')') {
                throw std::runtime_error("Mismatched parentheses");
            }
            ++pos;
            return result;
        }
        if (c == '-') {
            ++pos;
            return -parseFactor();
        }
        if (std::isdigit(c) || c == '.') {
            return parseNumber();
        }
        std::string token = parseToken();
        if (functions.count(token)) {
            skipWhitespace();
            if (pos >= expression.length() || expression[pos] != '(') {
                throw std::runtime_error("Expected '(' after function " + token);
            }
            ++pos;
            long double arg = parseExpression();
            skipWhitespace();
            if (pos >= expression.length() || expression[pos] != ')') {
                throw std::runtime_error("Expected ')' after function argument");
            }
            ++pos;
            long double result = functions[token](arg);
            if (std::isinf(result) || std::isnan(result)) {
                throw std::runtime_error("Numerical overflow in function " + token);
            }
            return result;
        }
        if (binary_functions.count(token)) {
            skipWhitespace();
            if (pos >= expression.length() || expression[pos] != '(') {
                throw std::runtime_error("Expected '(' after binary function " + token);
            }
            ++pos;
            long double arg1 = parseExpression();
            skipWhitespace();
            if (pos >= expression.length() || expression[pos] != ',') {
                throw std::runtime_error("Expected ',' in binary function " + token);
            }
            ++pos;
            long double arg2 = parseExpression();
            skipWhitespace();
            if (pos >= expression.length() || expression[pos] != ')') {
                throw std::runtime_error("Expected ')' after binary function arguments");
            }
            ++pos;
            long double result = binary_functions[token](arg1, arg2);
            if (std::isinf(result) || std::isnan(result)) {
                throw std::runtime_error("Numerical overflow in binary function " + token);
            }
            return result;
        }
        if (vars.count(token)) {
            return vars[token];
        }
        throw std::runtime_error("Unknown variable or function: " + token);
    }

    long double parseNumber() {
        skipWhitespace();
        size_t start = pos;
        while (pos < expression.length() && (std::isdigit(expression[pos]) || expression[pos] == '.' || expression[pos] == 'e' || expression[pos] == '+' || expression[pos] == '-')) {
            ++pos;
        }
        std::string num_str = expression.substr(start, pos - start);
        try {
            size_t idx;
            long double value = std::stold(num_str, &idx);
            if (idx != num_str.length()) {
                throw std::runtime_error("Invalid number format: " + num_str);
            }
            return value;
        } catch (...) {
            throw std::runtime_error("Invalid number: " + num_str);
        }
    }

    std::string parseToken() {
        skipWhitespace();
        size_t start = pos;
        while (pos < expression.length() && (std::isalnum(expression[pos]) || expression[pos] == '_')) {
            ++pos;
        }
        if (start == pos) {
            throw std::runtime_error("Expected token at position " + std::to_string(pos));
        }
        return expression.substr(start, pos - start);
    }

public:
    ExpressionParser() {
        functions["sin"] = [](long double x) { return std::sin(x); };
        functions["cos"] = [](long double x) { return std::cos(x); };
        functions["sqrt"] = [](long double x) { return x >= 0 ? std::sqrt(x) : throw std::runtime_error("sqrt of negative number"); };
        functions["exp"] = [](long double x) { return x <= 700.0L ? std::exp(x) : throw std::runtime_error("exp overflow"); };
        functions["log"] = [](long double x) { return x > 0 ? std::log(x) : throw std::runtime_error("log of non-positive number"); };
        functions["tan"] = [](long double x) { return std::tan(x); };
        functions["acos"] = [](long double x) { return x >= -1 && x <= 1 ? std::acos(x) : throw std::runtime_error("acos domain error"); };
        functions["asin"] = [](long double x) { return x >= -1 && x <= 1 ? std::asin(x) : throw std::runtime_error("asin domain error"); };
        functions["atan"] = [](long double x) { return std::atan(x); };
        binary_functions["max"] = [](long double x, long double y) { return std::max(x, y); };
        binary_functions["min"] = [](long double x, long double y) { return std::min(x, y); };
        binary_functions["atan2"] = [](long double y, long double x) { return std::atan2(y, x); };
    }

    long double evaluate(const std::string& expr, const std::map<std::string, long double>& vars) {
        this->vars = vars;
        pos = 0;
        expression = preprocess(expr) + " ";
        try {
            long double result = parseExpression();
            if (pos < expression.length() && !std::isspace(expression[pos])) {
                throw std::runtime_error("Unexpected character at position " + std::to_string(pos));
            }
            return result;
        } catch (const std::exception& e) {
            throw std::runtime_error("Evaluation error: " + std::string(e.what()));
        }
    }

    bool validate(const std::string& expr) const {
        for (char c : expr) {
            if (!std::isalnum(c) && c != '_' && c != '+' && c != '-' && c != '*' && c != '/' &&
                c != '^' && c != '(' && c != ')' && c != '.' && c != ' ' && c != '=' && c != ',' &&
                c != 'e') {
                return false;
            }
        }
        return true;
    }
};

// Simulation class
class Simulation {
private:
    friend class CustomEquation;
    friend class QuantumEquation;

public:
    Simulation(int grid_size = 64, long double delta_x = 0.1L, long double delta_t = 0.01L);
    void load_equations();
    void init(int grid_size, long double dx_val, long double dt_val, int num_particles);
    void run(int steps, const std::string& eq_name);
    std::string output(const std::string& eq_name) const;
    void setVisParameters(const std::string& eq_name, long double scale, long double color_intensity);
    void setVisScale(const std::string& eq_name, long double scale);
    void setVisColorIntensity(const std::string& eq_name, long double intensity);
    void reset();
    void step();
    void setup();
    void prev_equation();
    void next_equation();
    void compute();

    long double t;
    long double dt;
    int size;
    long double dx;
    std::vector<Particle> particles;
    std::vector<std::vector<std::vector<Tensor>>> ricci;
    std::vector<std::vector<std::vector<Tensor>>> divergence;
    std::vector<std::vector<std::vector<long double>>> computed_scalar;
    std::vector<std::vector<std::vector<std::vector<long double>>>> scalar_4d;
    std::vector<std::unique_ptr<Equation>> equations;
    size_t current_equation;
};

// Inline Simulation methods
inline Simulation::Simulation(int grid_size, long double delta_x, long double delta_t)
    : t(0.0L), dt(delta_t), size(grid_size), dx(delta_x), current_equation(0) {
    ricci.resize(size, std::vector<std::vector<Tensor>>(size, std::vector<Tensor>(size)));
    divergence.resize(size, std::vector<std::vector<Tensor>>(size, std::vector<Tensor>(size)));
    computed_scalar.resize(size, std::vector<std::vector<long double>>(size, std::vector<long double>(size, 0.0L)));
    scalar_4d.resize(size, std::vector<std::vector<std::vector<long double>>>(
        size, std::vector<std::vector<long double>>(size, std::vector<long double>(size, 0.0L))));
}

inline void Simulation::step() {
    extern Logger logger;
    if (!equations.empty() && current_equation < equations.size()) {
        equations[current_equation]->compute(*this);
        t += dt;
    }
}

inline void Simulation::setup() {
    reset();
    load_equations();
    if (!equations.empty() && current_equation < equations.size()) {
        equations[current_equation]->compute(*this);
    }
}

inline void Simulation::prev_equation() {
    if (!equations.empty()) {
        current_equation = (current_equation == 0) ? equations.size() - 1 : current_equation - 1;
        compute();
    }
}

inline void Simulation::next_equation() {
    if (!equations.empty()) {
        current_equation = (current_equation + 1) % equations.size();
        compute();
    }
}

inline void Simulation::compute() {
    if (!equations.empty() && current_equation < equations.size()) {
        equations[current_equation]->compute(*this);
    }
}

class CustomEquation : public Equation {
private:
    std::string eq_name;
    std::string equation;
    std::map<std::string, long double> parameters;
    std::map<std::string, std::string> derived_parameters;
    bool is_particle_based;
    bool compute_divergence; // Added to support TensorEquation functionality
    bool valid;
    ExpressionParser parser;
    long double scale, mass, c, G, mass1, mass2, r;
    long double wave_speed, freq, h, planck_freq, a, T, M, radius;
    long double l, m, n, l_quant, m_quant, gw_amp, gw_freq;
    long double E, B, kg_mass, spin, viscosity, temp, dim, alpha;
    long double H0, H1D, D_t, eps_5D, I0, D0, FiveD0, k, w, a0;
    long double vis_scale, vis_color_intensity;
    const long double max_ld = std::numeric_limits<long double>::max();
    const long double min_ld = -std::numeric_limits<long double>::max();
    const long double duration = 30.0L;
    std::chrono::steady_clock::time_point start_time;

    long double compute_sum_field(long double t, const std::map<std::string, long double>& vars, const std::string& sum_expr) const {
        long double field_value = 0.0L;
        long double progress = std::min(t / duration, 1.0L);
        long double range_factor = min_ld + progress * (max_ld - min_ld);

        long double H = vars.at("H");
        long double alpha = vars.count("alpha") ? vars.at("alpha") : 0.01L;
        long double eps = vars.at("eps");
        long double L0 = vars.at("L0");
        long double scale_factor = vars.count("scale") ? vars.at("scale") : 1.0L;
        int d_max = static_cast<int>(vars.count("d_max") ? vars.at("d_max") : 100.0L);
        int n_max = static_cast<int>(vars.count("n_max") ? vars.at("n_max") : 50.0L);
        long double c_squared = vars.count("c") ? vars.at("c") * vars.at("c") : 8.987551789e16L;
        const long double convergence_threshold = vars.count("conv_th") ? vars.at("conv_th") : 1e-20L;

        d_max += static_cast<int>(progress * 900);
        n_max += static_cast<int>(progress * 450);

        for (int d = 5; d <= d_max; ++d) {
            long double Hd_t = H * std::exp(-alpha * d * t);
            long double d_factor = 1.0L / (static_cast<long double>(d) * d * d * d);
            long double inner_sum = 0.0L;

            for (int n = 2; n <= n_max; ++n) {
                long double epsilon_n_d = eps / (n * n);
                long double L_n_d = L0 / std::sqrt(static_cast<long double>(d * n));
                long double term = epsilon_n_d * c_squared * L_n_d * L_n_d;
                if (!std::isinf(term) && !std::isnan(term) && std::abs(term) > convergence_threshold) {
                    inner_sum += term;
                } else {
                    break;
                }
            }

            long double dimension_term = Hd_t * inner_sum * d_factor;
            if (!std::isinf(dimension_term) && !std::isnan(dimension_term)) {
                field_value += dimension_term;
            } else {
                break;
            }
        }

        field_value *= scale_factor * range_factor;
        const long double max_vis_value = 1.0e10L;
        if (std::abs(field_value) > max_vis_value) {
            field_value = std::copysign(max_vis_value, field_value);
        }

        return field_value;
    }

public:
    CustomEquation(const std::string& name, const std::string& eq, const std::map<std::string, long double>& params,
                   const std::map<std::string, std::string>& derived_params, bool is_particle_based, bool compute_divergence = false)
        : eq_name(name), equation(eq), parameters(params), derived_parameters(derived_params),
          is_particle_based(is_particle_based), compute_divergence(compute_divergence), valid(parser.validate(eq)),
          scale(1.0L), mass(1.0L), c(299792458.0L), G(6.67430e-11L), mass1(1.0L), mass2(1.0L), r(1.0L),
          wave_speed(343.0L), freq(1.0L), h(6.62607015e-34L), planck_freq(1.0e15L), a(1.0L), T(1.0L),
          M(1.0L), radius(1.0L), l(0.0L), m(0.0L), n(1.0L), l_quant(0.0L), m_quant(0.0L),
          gw_amp(1.0L), gw_freq(1.0L), E(1.0L), B(1.0L), kg_mass(1.0L), spin(0.5L), viscosity(0.01L),
          temp(5000.0L), dim(10.0L), alpha(0.01L), H0(1.0e-18L), H1D(1.0L), D_t(0.1L), eps_5D(0.1L),
          I0(0.0L), D0(0.1L), FiveD0(1.0L), k(1.0L), w(1.0L), a0(5.29e-11L), vis_scale(1.5L), vis_color_intensity(1.4L),
          start_time(std::chrono::steady_clock::now()) {
        extern Logger logger;
        for (const auto& param : params) {
            const std::string& key = param.first;
            if (key == "scale") scale = param.second;
            else if (key == "mass") mass = param.second;
            else if (key == "c") c = param.second;
            else if (key == "G") G = param.second;
            else if (key == "mass1") mass1 = param.second;
            else if (key == "mass2") mass2 = param.second;
            else if (key == "r") r = param.second;
            else if (key == "wave_speed") wave_speed = param.second;
            else if (key == "freq") freq = param.second;
            else if (key == "h") h = param.second;
            else if (key == "planck_freq") planck_freq = param.second;
            else if (key == "a") a = param.second;
            else if (key == "T") T = param.second;
            else if (key == "M") M = param.second;
            else if (key == "radius") radius = param.second;
            else if (key == "l") l = param.second;
            else if (key == "m") m = param.second;
            else if (key == "n") n = param.second;
            else if (key == "l_quant") l_quant = param.second;
            else if (key == "m_quant") m_quant = param.second;
            else if (key == "amp") gw_amp = param.second;
            else if (key == "gw_freq") gw_freq = param.second;
            else if (key == "E") E = param.second;
            else if (key == "B") B = param.second;
            else if (key == "kg_mass") kg_mass = param.second;
            else if (key == "spin") spin = param.second;
            else if (key == "viscosity") viscosity = param.second;
            else if (key == "temp") temp = param.second;
            else if (key == "dim") dim = param.second;
            else if (key == "alpha") alpha = param.second;
            else if (key == "H0") H0 = param.second;
            else if (key == "H1D") H1D = param.second;
            else if (key == "D_t") D_t = param.second;
            else if (key == "eps_5D") eps_5D = param.second;
            else if (key == "I0") I0 = param.second;
            else if (key == "D0") D0 = param.second;
            else if (key == "FiveD0") FiveD0 = param.second;
            else if (key == "k") k = param.second;
            else if (key == "w") w = param.second;
            else if (key == "a0") a0 = param.second;
            else if (key == "vis_scale") vis_scale = param.second;
            else if (key == "vis_color_intensity") vis_color_intensity = param.second;
            else {
                parameters[key] = param.second;
                logger.log("INFO: Accepted parameter: " + key + "=" + std::to_string(param.second));
            }
        }
        if (!valid) {
            logger.log("Invalid equation syntax: " + eq);
        }
    }

    bool is_valid() const override;
    void compute(Simulation& sim) override;
    std::string output(const Simulation& sim) const override;
    std::string name() const override;
    void setVisScale(long double s) override;
    long double getVisScale() const override;
    void setVisColorIntensity(long double c) override;
    long double getVisColorIntensity() const override;
};

class QuantumEquation : public Equation {
private:
    std::string eq_name;
    long double hbar, m, a0;
    long double vis_scale, vis_color_intensity;

public:
    QuantumEquation(const std::string& name, const std::map<std::string, long double>& params)
        : eq_name(name), hbar(1.05457182e-34L), m(9.1093837e-31L), a0(5.29e-11L),
          vis_scale(1.5L), vis_color_intensity(1.4L) {
        extern Logger logger;
        for (const auto& param : params) {
            const std::string& key = param.first;
            if (key == "hbar") hbar = param.second;
            else if (key == "m") m = param.second;
            else if (key == "a0") a0 = param.second;
            else if (key == "vis_scale") vis_scale = param.second;
            else if (key == "vis_color_intensity") vis_color_intensity = param.second;
            else {
                logger.log("INFO: Accepted parameter: " + key + "=" + std::to_string(param.second));
            }
        }
    }

    bool is_valid() const override;
    void compute(Simulation& sim) override;
    std::string output(const Simulation& sim) const override;
    std::string name() const override;
    void setVisScale(long double s) override;
    long double getVisScale() const override;
    void setVisColorIntensity(long double c) override;
    long double getVisColorIntensity() const override;
};

#endif