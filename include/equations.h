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
    std::vector<std::string> buffer; // Corrected syntax
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
    friend class TensorEquation;
    friend class QuantumEquation;
    friend class OriginalEquation;
    friend class RefinedEquation;

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
                   const std::map<std::string, std::string>& derived_params, bool is_particle_based)
        : eq_name(name), equation(eq), parameters(params), derived_parameters(derived_params),
          is_particle_based(is_particle_based), valid(parser.validate(eq)),
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

    bool is_valid() const override {
        return valid;
    }

    void compute(Simulation& sim) override {
        extern Logger logger;
        std::map<std::string, long double> vars = {
            {"H", H1D}, {"eps", eps_5D}, {"L0", 1.0L}, {"scale", scale}, {"alpha", alpha}
        };
        #pragma omp parallel for collapse(3)
        for (int i = 0; i < sim.size; ++i) {
            for (int j = 0; j < sim.size; ++j) {
                for (int k = 0; k < sim.size; ++k) {
                    long double x = static_cast<long double>(i - sim.size / 2) * sim.dx;
                    long double y = static_cast<long double>(j - sim.size / 2) * sim.dx;
                    long double z = static_cast<long double>(k - sim.size / 2) * sim.dx;
                    long double r = std::sqrt(x * x + y * y + z * z) + 1e-6L;
                    long double scalar = compute_sum_field(sim.t, vars, equation);
                    scalar = std::tanh(scalar / vis_color_intensity) * vis_scale;
                    #pragma omp critical
                    {
                        sim.computed_scalar[i][j][k] = scalar;
                        sim.ricci[i][j][k].set_diagonal(scalar);
                        for (int w = 0; w < sim.size; ++w) {
                            long double w_coord = static_cast<long double>(w - sim.size / 2) * sim.dx;
                            sim.scalar_4d[i][j][k][w] = scalar * (1.0L + 0.1L * std::cos(w_coord + sim.t));
                        }
                    }
                }
            }
        }
        sim.t += sim.dt;
        logger.flush();
    }

    std::string output(const Simulation& sim) const override {
        extern Logger logger;
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(6);
        oss << "# Equation: " << name() << "\n";
        oss << "# Time: " << sim.t << "\n";
        if (is_particle_based) {
            oss << "# Particles (index,x,y,z,vx,vy,vz,ax,ay,az)\n";
            long double total_Ke = 0.0L;
            for (size_t i = 0; i < sim.particles.size(); ++i) {
                const auto& p = sim.particles[i];
                oss << i << "," << p.x << "," << p.y << "," << p.z << ","
                    << p.vx << "," << p.vy << "," << p.vz << ","
                    << p.ax << "," << p.ay << "," << p.az << "\n";
                long double speed = std::sqrt(p.vx * p.vx + p.vy * p.vy + p.vz * p.vz);
                total_Ke += 0.5L * mass * speed * speed;
            }
            oss << "# Total Kinetic Energy: " << total_Ke << "\n";
        } else {
            oss << "# Scalar Field and Gradients (i,j,k,x,y,z,scalar,grad_x,grad_y,grad_z)\n";
            long double sum_scalar = 0.0L;
            int count = 0;
            for (int i = 0; i < sim.size; ++i) {
                for (int j = 0; j < sim.size; ++j) {
                    for (int k = 0; k < sim.size; ++k) {
                        long double x = static_cast<long double>(i - sim.size / 2) * sim.dx;
                        long double y = static_cast<long double>(j - sim.size / 2) * sim.dx;
                        long double z = static_cast<long double>(k - sim.size / 2) * sim.dx;
                        long double scalar = sim.computed_scalar[i][j][k];
                        sum_scalar += scalar;
                        count++;
                        long double grad_x = (i > 0 && i < sim.size - 1) ? (sim.computed_scalar[i + 1][j][k] - sim.computed_scalar[i - 1][j][k]) / (2 * sim.dx) : 0.0L;
                        long double grad_y = (j > 0 && j < sim.size - 1) ? (sim.computed_scalar[i][j + 1][k] - sim.computed_scalar[i][j - 1][k]) / (2 * sim.dx) : 0.0L;
                        long double grad_z = (k > 0 && k < sim.size - 1) ? (sim.computed_scalar[i][j][k + 1] - sim.computed_scalar[i][j][k - 1]) / (2 * sim.dx) : 0.0L;
                        oss << i << "," << j << "," << k << "," << x << "," << y << "," << z << ","
                            << scalar << "," << grad_x << "," << grad_y << "," << grad_z << "\n";
                    }
                }
            }
            oss << "# Average Scalar: " << (count > 0 ? sum_scalar / count : 0.0L) << "\n";
            oss << "# Total Scalar (Integral): " << (sum_scalar * sim.dx * sim.dx * sim.dx) << "\n";
        }
        return oss.str();
    }

    std::string name() const override {
        return eq_name;
    }

    void setVisScale(long double s) override {
        vis_scale = std::max(0.1L, std::min(10.0L, s));
    }

    long double getVisScale() const override {
        return vis_scale;
    }

    void setVisColorIntensity(long double c) override {
        vis_color_intensity = std::max(0.1L, std::min(2.0L, c));
    }

    long double getVisColorIntensity() const override {
        return vis_color_intensity;
    }
};

class TensorEquation : public Equation {
private:
    std::string eq_name;
    std::string equation;
    std::map<std::string, long double> parameters;
    std::map<std::string, std::string> derived_parameters;
    long double H1D, D_t_coeff, eps_5D, alpha, I0, D0, FiveD0;
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
    TensorEquation(const std::string& name, const std::string& eq, const std::map<std::string, long double>& params,
                   const std::map<std::string, std::string>& derived_params)
        : eq_name(name), equation(eq), parameters(params), derived_parameters(derived_params),
          H1D(1.0L), D_t_coeff(0.1L), eps_5D(0.1L), alpha(0.01L), I0(0.0L), D0(0.1L), FiveD0(1.0L),
          vis_scale(1.5L), vis_color_intensity(1.4L), start_time(std::chrono::steady_clock::now()) {
        extern Logger logger;
        for (const auto& param : params) {
            const std::string& key = param.first;
            if (key == "H1D") H1D = param.second;
            else if (key == "D_t") D_t_coeff = param.second;
            else if (key == "eps_5D") eps_5D = param.second;
            else if (key == "alpha") alpha = param.second;
            else if (key == "I0") I0 = param.second;
            else if (key == "D0") D0 = param.second;
            else if (key == "FiveD0") FiveD0 = param.second;
            else if (key == "vis_scale") vis_scale = param.second;
            else if (key == "vis_color_intensity") vis_color_intensity = param.second;
            else {
                logger.log("INFO: Accepted parameter: " + key + "=" + std::to_string(param.second));
            }
        }
    }

    bool is_valid() const override {
        ExpressionParser parser;
        return parser.validate(equation);
    }

    void compute(Simulation& sim) override {
        extern Logger logger;
        std::map<std::string, long double> vars = {
            {"H", H1D}, {"eps", eps_5D}, {"L0", FiveD0}, {"scale", vis_scale}, {"alpha", alpha}
        };
        #pragma omp parallel for collapse(3)
        for (int i = 0; i < sim.size; ++i) {
            for (int j = 0; j < sim.size; ++j) {
                for (int k = 0; k < sim.size; ++k) {
                    long double x = static_cast<long double>(i - sim.size / 2) * sim.dx;
                    long double y = static_cast<long double>(j - sim.size / 2) * sim.dx;
                    long double z = static_cast<long double>(k - sim.size / 2) * sim.dx;
                    long double r = std::sqrt(x * x + y * y + z * z) + 1e-6L;
                    long double scalar = compute_sum_field(sim.t, vars, equation);
                    scalar = std::tanh(scalar / vis_color_intensity) * vis_scale;
                    #pragma omp critical
                    {
                        sim.computed_scalar[i][j][k] = scalar;
                        sim.ricci[i][j][k].set_diagonal(scalar + I0);
                        sim.divergence[i][j][k].set_diagonal(D0 * std::sin(sim.t));
                        for (int w = 0; w < sim.size; ++w) {
                            long double w_coord = static_cast<long double>(w - sim.size / 2) * sim.dx;
                            sim.scalar_4d[i][j][k][w] = scalar * (1.0L + D_t_coeff * std::cos(w_coord + sim.t));
                        }
                    }
                }
            }
        }
        sim.t += sim.dt;
        logger.flush();
    }

    std::string output(const Simulation& sim) const override {
        extern Logger logger;
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(6);
        oss << "# Equation: " << name() << "\n";
        oss << "# Time: " << sim.t << "\n";
        oss << "# Scalar Field and Gradients (i,j,k,x,y,z,scalar,grad_x,grad_y,grad_z)\n";
        long double sum_scalar = 0.0L;
        int count = 0;
        for (int i = 0; i < sim.size; ++i) {
            for (int j = 0; j < sim.size; ++j) {
                for (int k = 0; k < sim.size; ++k) {
                    long double x = static_cast<long double>(i - sim.size / 2) * sim.dx;
                    long double y = static_cast<long double>(j - sim.size / 2) * sim.dx;
                    long double z = static_cast<long double>(k - sim.size / 2) * sim.dx;
                    long double scalar = sim.computed_scalar[i][j][k];
                    sum_scalar += scalar;
                    count++;
                    long double grad_x = (i > 0 && i < sim.size - 1) ? (sim.computed_scalar[i + 1][j][k] - sim.computed_scalar[i - 1][j][k]) / (2 * sim.dx) : 0.0L;
                    long double grad_y = (j > 0 && j < sim.size - 1) ? (sim.computed_scalar[i][j + 1][k] - sim.computed_scalar[i][j - 1][k]) / (2 * sim.dx) : 0.0L;
                    long double grad_z = (k > 0 && k < sim.size - 1) ? (sim.computed_scalar[i][j][k + 1] - sim.computed_scalar[i][j][k - 1]) / (2 * sim.dx) : 0.0L;
                    oss << i << "," << j << "," << k << "," << x << "," << y << "," << z << ","
                        << scalar << "," << grad_x << "," << grad_y << "," << grad_z << "\n";
                }
            }
        }
        oss << "# Average Scalar: " << (count > 0 ? sum_scalar / count : 0.0L) << "\n";
        oss << "# Total Scalar (Integral): " << (sum_scalar * sim.dx * sim.dx * sim.dx) << "\n";
        return oss.str();
    }

    std::string name() const override {
        return eq_name;
    }

    void setVisScale(long double s) override {
        vis_scale = std::max(0.1L, std::min(10.0L, s));
    }

    long double getVisScale() const override {
        return vis_scale;
    }

    void setVisColorIntensity(long double c) override {
        vis_color_intensity = std::max(0.1L, std::min(2.0L, c));
    }

    long double getVisColorIntensity() const override {
        return vis_color_intensity;
    }
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

    bool is_valid() const override {
        return true;
    }

    void compute(Simulation& sim) override {
        extern Logger logger;
        #pragma omp parallel for collapse(3)
        for (int i = 0; i < sim.size; ++i) {
            for (int j = 0; j < sim.size; ++j) {
                for (int k = 0; k < sim.size; ++k) {
                    long double x = static_cast<long double>(i - sim.size / 2) * sim.dx;
                    long double y = static_cast<long double>(j - sim.size / 2) * sim.dx;
                    long double z = static_cast<long double>(k - sim.size / 2) * sim.dx;
                    long double r = std::sqrt(x * x + y * y + z * z) + 1e-6L;
                    long double psi = std::exp(-r / a0) / (std::sqrt(M_PI) * a0 * a0);
                    long double psi_squared = psi * psi;
                    long double vis_value = psi_squared * vis_scale;
                    vis_value = std::tanh(vis_value / vis_color_intensity);
                    #pragma omp critical
                    {
                        sim.computed_scalar[i][j][k] = vis_value;
                        sim.ricci[i][j][k].set_diagonal(vis_value);
                        for (int w = 0; w < sim.size; ++w) {
                            long double w_coord = static_cast<long double>(w - sim.size / 2) * sim.dx;
                            sim.scalar_4d[i][j][k][w] = vis_value * (1.0L + 0.1L * std::cos(w_coord + sim.t));
                        }
                    }
                }
            }
        }
        sim.t += sim.dt;
        logger.flush();
    }

    std::string output(const Simulation& sim) const override {
        extern Logger logger;
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(6);
        oss << "# Equation: " << name() << "\n";
        oss << "# Time: " << sim.t << "\n";
        oss << "# Scalar Field (Probability Density) and Gradients (i,j,k,x,y,z,scalar,grad_x,grad_y,grad_z)\n";
        long double sum_scalar = 0.0L;
        int count = 0;
        for (int i = 0; i < sim.size; ++i) {
            for (int j = 0; j < sim.size; ++j) {
                for (int k = 0; k < sim.size; ++k) {
                    long double x = static_cast<long double>(i - sim.size / 2) * sim.dx;
                    long double y = static_cast<long double>(j - sim.size / 2) * sim.dx;
                    long double z = static_cast<long double>(k - sim.size / 2) * sim.dx;
                    long double scalar = sim.computed_scalar[i][j][k];
                    sum_scalar += scalar;
                    count++;
                    long double grad_x = (i > 0 && i < sim.size - 1) ? (sim.computed_scalar[i + 1][j][k] - sim.computed_scalar[i - 1][j][k]) / (2 * sim.dx) : 0.0L;
                    long double grad_y = (j > 0 && j < sim.size - 1) ? (sim.computed_scalar[i][j + 1][k] - sim.computed_scalar[i][j - 1][k]) / (2 * sim.dx) : 0.0L;
                    long double grad_z = (k > 0 && k < sim.size - 1) ? (sim.computed_scalar[i][j][k + 1] - sim.computed_scalar[i][j][k - 1]) / (2 * sim.dx) : 0.0L;
                    oss << i << "," << j << "," << k << "," << x << "," << y << "," << z << ","
                        << scalar << "," << grad_x << "," << grad_y << "," << grad_z << "\n";
                }
            }
        }
        oss << "# Average Scalar: " << (count > 0 ? sum_scalar / count : 0.0L) << "\n";
        oss << "# Total Scalar (Integral): " << (sum_scalar * sim.dx * sim.dx * sim.dx) << "\n";
        return oss.str();
    }

    std::string name() const override {
        return eq_name;
    }

    void setVisScale(long double s) override {
        vis_scale = std::max(0.1L, std::min(10.0L, s));
    }

    long double getVisScale() const override {
        return vis_scale;
    }

    void setVisColorIntensity(long double c) override {
        vis_color_intensity = std::max(0.1L, std::min(2.0L, c));
    }

    long double getVisColorIntensity() const override {
        return vis_color_intensity;
    }
};

class OriginalEquation : public Equation {
private:
    std::string eq_name;
    long double H1D, D_t_coeff, eps_5D;
    long double vis_scale, vis_color_intensity;

public:
    OriginalEquation(const std::string& name, const std::map<std::string, long double>& params)
        : eq_name(name), H1D(1.0L), D_t_coeff(0.1L), eps_5D(0.1L), vis_scale(1.5L), vis_color_intensity(1.4L) {
        extern Logger logger;
        for (const auto& param : params) {
            const std::string& key = param.first;
            if (key == "H1D") H1D = param.second;
            else if (key == "D_t") D_t_coeff = param.second;
            else if (key == "eps_5D") eps_5D = param.second;
            else if (key == "vis_scale") vis_scale = param.second;
            else if (key == "vis_color_intensity") vis_color_intensity = param.second;
            else {
                logger.log("INFO: Accepted parameter: " + key + "=" + std::to_string(param.second));
            }
        }
    }

    bool is_valid() const override {
        return true;
    }

    void compute(Simulation& sim) override {
        extern Logger logger;
        long double D_t = D_t_coeff * std::sin(sim.t);
        long double eps = eps_5D * std::cos(sim.t);

        #pragma omp parallel for collapse(3)
        for (int i = 1; i < sim.size - 1; ++i) {
            for (int j = 1; j < sim.size - 1; ++j) {
                for (int k = 1; k < sim.size - 1; ++k) {
                    Tensor div;
                    for (int nu = 0; nu < 3; ++nu) {
                        for (int mu = 0; mu < 3; ++mu) {
                            long double dR = 0.0L;
                            if (mu == 0) {
                                dR = (sim.ricci[i + 1][j][k].data[mu][nu] - sim.ricci[i - 1][j][k].data[mu][nu]) / (2 * sim.dx);
                            } else if (mu == 1) {
                                dR = (sim.ricci[i][j + 1][k].data[mu][nu] - sim.ricci[i][j - 1][k].data[mu][nu]) / (2 * sim.dx);
                            } else if (mu == 2) {
                                dR = (sim.ricci[i][j][k + 1].data[mu][nu] - sim.ricci[i][j][k - 1].data[mu][nu]) / (2 * sim.dx);
                            }
                            div.data[mu][nu] = dR;
                        }
                    }

                    Tensor rhs;
                    rhs.set_diagonal(H1D);
                    for (int m = 0; m < 3; ++m) {
                        rhs.data[m][m] += sim.ricci[i][j][k].data[m][m];
                        rhs.data[m][m] -= 4 * D_t;
                        rhs.data[m][m] += eps;
                    }

                    long double scalar = 0.0L;
                    for (int m = 0; m < 3; ++m) {
                        scalar += div.data[m][m];
                    }
                    scalar = std::sqrt(std::abs(scalar));
                    long double vis_value = scalar * vis_scale;
                    vis_value = std::tanh(vis_value / vis_color_intensity);

                    #pragma omp critical
                    {
                        sim.divergence[i][j][k] = div;
                        sim.computed_scalar[i][j][k] = vis_value;
                        for (int w = 0; w < sim.size; ++w) {
                            long double w_coord = static_cast<long double>(w - sim.size / 2) * sim.dx;
                            sim.scalar_4d[i][j][k][w] = vis_value * (1.0L + 0.1L * std::sin(w_coord + sim.t));
                        }
                    }
                }
            }
        }
        sim.t += sim.dt;
        logger.flush();
    }

    std::string output(const Simulation& sim) const override {
        extern Logger logger;
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(6);
        oss << "# Equation: " << name() << "\n";
        oss << "# Time: " << sim.t << "\n";
        oss << "# Scalar Field and Gradients (i,j,k,x,y,z,scalar,grad_x,grad_y,grad_z)\n";
        long double sum_scalar = 0.0L;
        int count = 0;
        for (int i = 0; i < sim.size; ++i) {
            for (int j = 0; j < sim.size; ++j) {
                for (int k = 0; k < sim.size; ++k) {
                    long double x = static_cast<long double>(i - sim.size / 2) * sim.dx;
                    long double y = static_cast<long double>(j - sim.size / 2) * sim.dx;
                    long double z = static_cast<long double>(k - sim.size / 2) * sim.dx;
                    long double scalar = sim.computed_scalar[i][j][k];
                    sum_scalar += scalar;
                    count++;
                    long double grad_x = (i > 0 && i < sim.size - 1) ? (sim.computed_scalar[i + 1][j][k] - sim.computed_scalar[i - 1][j][k]) / (2 * sim.dx) : 0.0L;
                    long double grad_y = (j > 0 && j < sim.size - 1) ? (sim.computed_scalar[i][j + 1][k] - sim.computed_scalar[i][j - 1][k]) / (2 * sim.dx) : 0.0L;
                    long double grad_z = (k > 0 && k < sim.size - 1) ? (sim.computed_scalar[i][j][k + 1] - sim.computed_scalar[i][j][k - 1]) / (2 * sim.dx) : 0.0L;
                    oss << i << "," << j << "," << k << "," << x << "," << y << "," << z << ","
                        << scalar << "," << grad_x << "," << grad_y << "," << grad_z << "\n";
                }
            }
        }
        oss << "# Average Scalar: " << (count > 0 ? sum_scalar / count : 0.0L) << "\n";
        oss << "# Total Scalar (Integral): " << (sum_scalar * sim.dx * sim.dx * sim.dx) << "\n";
        return oss.str();
    }

    std::string name() const override {
        return eq_name;
    }

    void setVisScale(long double s) override {
        vis_scale = std::max(0.1L, std::min(10.0L, s));
    }

    long double getVisScale() const override {
        return vis_scale;
    }

    void setVisColorIntensity(long double c) override {
        vis_color_intensity = std::max(0.1L, std::min(2.0L, c));
    }

    long double getVisColorIntensity() const override {
        return vis_color_intensity;
    }
};

class RefinedEquation : public Equation {
private:
    std::string eq_name;
    long double phi_base, G, k4, k5, eps_5D;
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
    RefinedEquation(const std::string& name, const std::map<std::string, long double>& params)
        : eq_name(name), phi_base(1.0L), G(1.0L), k4(0.1L), k5(0.01L), vis_scale(1.5L), vis_color_intensity(1.4L),
          start_time(std::chrono::steady_clock::now()) {
        extern Logger logger;
        for (const auto& param : params) {
            const std::string& key = param.first;
            if (key == "phi") phi_base = param.second;
            else if (key == "G") G = param.second;
            else if (key == "k4") k4 = param.second;
            else if (key == "k5") k5 = param.second;
            else if (key == "vis_scale") vis_scale = param.second;
            else if (key == "vis_color_intensity") vis_color_intensity = param.second;
            else {
                logger.log("INFO: Accepted parameter: " + key + "=" + std::to_string(param.second));
            }
        }
    }

    bool is_valid() const override {
        return true;
    }

    void compute(Simulation& sim) override {
        extern Logger logger;
        long double phi = phi_base * (1.0L + 0.1L * std::sin(sim.t));
        long double G_eff = G * (1.0L + k4 * std::cos(sim.t));
        long double k5_eff = k5 * std::exp(-0.01L * sim.t);

        #pragma omp parallel for collapse(3)
        for (int i = 1; i < sim.size - 1; ++i) {
            for (int j = 1; j < sim.size - 1; ++j) {
                for (int k = 1; k < sim.size - 1; ++k) {
                    long double x = static_cast<long double>(i - sim.size / 2) * sim.dx;
                    long double y = static_cast<long double>(j - sim.size / 2) * sim.dx;
                    long double z = static_cast<long double>(k - sim.size / 2) * sim.dx;
                    long double r = std::sqrt(x * x + y * y + z * z) + 1e-6L;

                    Tensor R;
                    for (int mu = 0; mu < 3; ++mu) {
                        for (int nu = 0; nu < 3; ++nu) {
                            long double dRdx = (sim.ricci[i + 1][j][k].data[mu][nu] - sim.ricci[i - 1][j][k].data[mu][nu]) / (2 * sim.dx);
                            long double dRdy = (sim.ricci[i][j + 1][k].data[mu][nu] - sim.ricci[i][j - 1][k].data[mu][nu]) / (2 * sim.dx);
                            long double dRdz = (sim.ricci[i][j][k + 1].data[mu][nu] - sim.ricci[i][j][k - 1].data[mu][nu]) / (2 * sim.dx);
                            R.data[mu][nu] = dRdx + dRdy + dRdz;
                        }
                    }

                    Tensor T; T.set_diagonal(0.0L);
                    Tensor F; F.set_diagonal(0.1L * std::sin(sim.t));
                    std::map<std::string, long double> vars = {{"H", phi}, {"eps", eps_5D}, {"L0", 1.0L}, {"scale", vis_scale}, {"alpha", 0.01L}};
                    long double eps_5D_val = compute_sum_field(sim.t, vars, "");
                    Tensor eps_5D; eps_5D.set_diagonal(eps_5D_val * k5_eff);

                    long double scalar = 0.0L;
                    for (int m = 0; m < 3; ++m) {
                        long double R_val = R.data[m][m];
                        long double T_val = T.data[m][m];
                        long double F_val = F.data[m][m];
                        long double eps_val = eps_5D.data[m][m];
                        scalar += (R_val - (phi * G_eff + T_val + F_val + eps_val)) * vis_scale;
                    }
                    scalar = std::tanh(std::abs(scalar) / vis_color_intensity);

                    #pragma omp critical
                    {
                        sim.ricci[i][j][k] = R;
                        sim.computed_scalar[i][j][k] = scalar;
                        for (int w = 0; w < sim.size; ++w) {
                            long double w_coord = static_cast<long double>(w - sim.size / 2) * sim.dx;
                            sim.scalar_4d[i][j][k][w] = scalar * (1.0L + 0.1L * std::cos(w_coord + sim.t));
                        }
                    }
                }
            }
        }
        sim.t += sim.dt;
        logger.flush();
    }

    std::string output(const Simulation& sim) const override {
        extern Logger logger;
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(6);
        oss << "# Equation: " << name() << "\n";
        oss << "# Time: " << sim.t << "\n";
        oss << "# Scalar Field and Gradients (i,j,k,x,y,z,scalar,grad_x,grad_y,grad_z)\n";
        long double sum_scalar = 0.0L;
        int count = 0;
        for (int i = 0; i < sim.size; ++i) {
            for (int j = 0; j < sim.size; ++j) {
                for (int k = 0; k < sim.size; ++k) {
                    long double x = static_cast<long double>(i - sim.size / 2) * sim.dx;
                    long double y = static_cast<long double>(j - sim.size / 2) * sim.dx;
                    long double z = static_cast<long double>(k - sim.size / 2) * sim.dx;
                    long double scalar = sim.computed_scalar[i][j][k];
                    sum_scalar += scalar;
                    count++;
                    long double grad_x = (i > 0 && i < sim.size - 1) ? 
                        (sim.computed_scalar[i + 1][j][k] - sim.computed_scalar[i - 1][j][k]) / (2 * sim.dx) : 0.0L;
                    long double grad_y = (j > 0 && j < sim.size - 1) ? 
                        (sim.computed_scalar[i][j + 1][k] - sim.computed_scalar[i][j - 1][k]) / (2 * sim.dx) : 0.0L;
                    long double grad_z = (k > 0 && k < sim.size - 1) ? 
                        (sim.computed_scalar[i][j][k + 1] - sim.computed_scalar[i][j][k - 1]) / (2 * sim.dx) : 0.0L;
                    oss << i << "," << j << "," << k << "," << x << "," << y << "," << z << ","
                        << scalar << "," << grad_x << "," << grad_y << "," << grad_z << "\n";
                }
            }
        }
        oss << "# Average Scalar: " << (count > 0 ? sum_scalar / count : 0.0L) << "\n";
        oss << "# Total Scalar (Integral): " << (sum_scalar * sim.dx * sim.dx * sim.dx) << "\n";
        return oss.str();
    }

    std::string name() const override {
        return eq_name;
    }

    void setVisScale(long double s) override {
        vis_scale = std::max(0.1L, std::min(10.0L, s));
    }

    long double getVisScale() const override {
        return vis_scale;
    }

    void setVisColorIntensity(long double c) override {
        vis_color_intensity = std::max(0.1L, std::min(2.0L, c));
    }

    long double getVisColorIntensity() const override {
        return vis_color_intensity;
    }
};

#endif