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
#include <omp.h>

// Declare log_message at the top
void log_message(const std::string& message);

class Tensor {
public:
    Tensor();
    void set_diagonal(double value);
    double data[3][3];
};

class Particle {
public:
    Particle();
    double x, y, z;
    double vx, vy, vz;
    double ax, ay, az;
};

class Simulation;

class Equation {
public:
    virtual ~Equation() {}
    virtual void compute(Simulation& sim) = 0;
    virtual std::string output(const Simulation& sim) const = 0;
    virtual std::string name() const = 0;
    virtual bool is_valid() const { return true; }
};

class Simulation {
public:
    Simulation();
    void initialize(int num_particles, int grid_size, double grid_dx, double dt_val);
    void load_equations(const std::string& filename);
    void compute();
    void save_output();
    void switch_equation(int index);
    void setup();
    void step();
    void prev_equation();
    void next_equation();

    double t;
    double dt;
    int current_equation;
    int size;
    double dx;
    std::vector<Particle> particles;
    std::vector<std::vector<std::vector<Tensor>>> ricci;
    std::vector<std::vector<std::vector<Tensor>>> divergence;
    std::vector<std::vector<std::vector<double>>> scalar;
    std::vector<std::vector<std::vector<double>>> computed_scalar;
    std::vector<std::vector<std::vector<std::vector<double>>>> scalar_4d;
    std::vector<Equation*> equations;
    std::set<std::string> output_written;
};

class ExpressionParser {
public:
    ExpressionParser() {
        // Unary functions
        functions["sin"] = [](double x) { return std::sin(x); };
        functions["cos"] = [](double x) { return std::cos(x); };
        functions["sqrt"] = [](double x) { return std::sqrt(x); };
        functions["exp"] = [](double x) { return std::exp(x); };
        functions["log"] = [](double x) { return std::log(x); };
        functions["tan"] = [](double x) { return std::tan(x); };
        functions["acos"] = [](double x) { return std::acos(x); }; // For theta in Hydrogen_Wavefunction

        // Binary functions (already defined in binary_functions map)
        binary_functions["max"] = [](double x, double y) { return std::max(x, y); };
        binary_functions["atan2"] = [](double y, double x) { return std::atan2(y, x); };
    }

    double evaluate(const std::string& expr, const std::map<std::string, double>& vars) {
        this->vars = vars;
        pos = 0;
        expression = preprocess(expr) + " ";
        try {
            double result = parseExpression();
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
                c != '^' && c != '(' && c != ')' && c != '.' && c != ' ' && c != '=' && c != ',') {
                return false;
            }
        }
        return true;
    }

private:
    std::map<std::string, std::function<double(double)>> functions;
    std::map<std::string, std::function<double(double, double)>> binary_functions;
    std::map<std::string, double> vars;
    std::string expression;
    size_t pos;

    std::string preprocess(const std::string& expr) {
        std::string result;
        bool last_was_digit = false;
        bool last_was_var = false;
        for (size_t i = 0; i < expr.length(); ++i) {
            char c = expr[i];
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
            last_was_digit = std::isdigit(c) || c == '.';
            last_was_var = std::isalpha(c);
        }
        return result;
    }

    void skipWhitespace() {
        while (pos < expression.length() && std::isspace(expression[pos])) {
            ++pos;
        }
    }

    double parseExpression() {
        double result = parseTerm();
        while (pos < expression.length()) {
            skipWhitespace();
            char op = expression[pos];
            if (op != '+' && op != '-') break;
            ++pos;
            double term = parseTerm();
            if (op == '+') result += term;
            else result -= term;
        }
        return result;
    }

    double parseTerm() {
        double result = parseFactor();
        while (pos < expression.length()) {
            skipWhitespace();
            char op = expression[pos];
            if (op == '*' && expression[pos + 1] == '*') {
                pos += 2;
                double factor = parseFactor();
                result = std::pow(result, factor);
            } else if (op != '*' && op != '/') {
                break;
            } else {
                ++pos;
                double factor = parseFactor();
                if (op == '*') result *= factor;
                else if (factor == 0) throw std::runtime_error("Division by zero");
                else result /= factor;
            }
        }
        return result;
    }

    double parseFactor() {
        skipWhitespace();
        if (pos >= expression.length()) throw std::runtime_error("Unexpected end of expression");

        char c = expression[pos];
        if (c == '(') {
            ++pos;
            double result = parseExpression();
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
            double arg = parseExpression();
            skipWhitespace();
            if (pos >= expression.length() || expression[pos] != ')') {
                throw std::runtime_error("Expected ')' after function argument");
            }
            ++pos;
            return functions[token](arg);
        }
        if (binary_functions.count(token)) {
            skipWhitespace();
            if (pos >= expression.length() || expression[pos] != '(') {
                throw std::runtime_error("Expected '(' after binary function " + token);
            }
            ++pos;
            double arg1 = parseExpression();
            skipWhitespace();
            if (pos >= expression.length() || expression[pos] != ',') {
                throw std::runtime_error("Expected ',' in binary function " + token);
            }
            ++pos;
            double arg2 = parseExpression();
            skipWhitespace();
            if (pos >= expression.length() || expression[pos] != ')') {
                throw std::runtime_error("Expected ')' after binary function arguments");
            }
            ++pos;
            return binary_functions[token](arg1, arg2);
        }
        if (vars.count(token)) {
            return vars[token];
        }
        throw std::runtime_error("Unknown variable or function: " + token);
    }

    double parseNumber() {
        skipWhitespace();
        size_t start = pos;
        while (pos < expression.length() && (std::isdigit(expression[pos]) || expression[pos] == '.' || expression[pos] == 'e' || expression[pos] == '+' || expression[pos] == '-')) {
            ++pos;
        }
        std::string num_str = expression.substr(start, pos - start);
        try {
            return std::stod(num_str);
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
        if (start == pos) throw std::runtime_error("Expected token at position " + std::to_string(pos));
        return expression.substr(start, pos - start);
    }
};

class CustomEquation : public Equation {
public:
    CustomEquation(const std::string& name, const std::string& eq, const std::map<std::string, double>& params, bool is_particle_based)
        : eq_name(name), equation(eq), parameters(params), is_particle_based(is_particle_based),
          scale(1.0), mass(1.0), c(299792458.0), G(6.67430e-11), mass1(1.0), mass2(1.0), r(1.0),
          wave_speed(343.0), freq(1.0), h(6.62607015e-34), planck_freq(1.0e15), a(1.0), T(1.0),
          M(1.0), radius(1.0), l(0.0), m(0.0), n(1.0), l_quant(0.0), m_quant(0.0),
          gw_amp(1.0), gw_freq(1.0), E(1.0), B(1.0), kg_mass(1.0), spin(0.5), viscosity(0.01),
          temp(5000.0), dim(10.0), alpha(1.0), H0(70.0), vis_scale(1.0), vis_color_intensity(1.0),
          valid(parser.validate(eq)) {
        // Initialize parameters from the provided map
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
            else if (key == "vis_scale") vis_scale = param.second;
            else if (key == "vis_color_intensity") vis_color_intensity = param.second;
            else {
                log_message("Warning: Unrecognized parameter '" + key + "' in equation '" + name + "'");
            }
        }
        if (!valid) {
            log_message("Invalid equation syntax: " + eq);
        }
    }

    bool is_valid() const override { return valid; }

    void compute(Simulation& sim) override {
        if (!valid) {
            sim.t += sim.dt;
            return;
        }
        if (is_particle_based) {
            for (auto& p : sim.particles) {
                std::map<std::string, double> vars = parameters;
                vars["x"] = p.x;
                vars["y"] = p.y;
                vars["z"] = p.z;
                vars["vx"] = p.vx;
                vars["vy"] = p.vy;
                vars["vz"] = p.vz;
                vars["ax"] = p.ax;
                vars["ay"] = p.ay;
                vars["az"] = p.az;
                vars["t"] = sim.t;
                vars["dt"] = sim.dt;
                vars["r"] = std::sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
                vars["theta"] = std::acos(p.z / std::max(vars["r"], 1e-6));
                vars["phi"] = std::atan2(p.y, p.x);
                // Add all member variables to vars
                vars["scale"] = scale;
                vars["mass"] = mass;
                vars["c"] = c;
                vars["G"] = G;
                vars["mass1"] = mass1;
                vars["mass2"] = mass2;
                vars["wave_speed"] = wave_speed;
                vars["freq"] = freq;
                vars["h"] = h;
                vars["planck_freq"] = planck_freq;
                vars["a"] = a;
                vars["T"] = T;
                vars["M"] = M;
                vars["radius"] = radius;
                vars["l"] = l;
                vars["m"] = m;
                vars["n"] = n;
                vars["l_quant"] = l_quant;
                vars["m_quant"] = m_quant;
                vars["gw_amp"] = gw_amp;
                vars["gw_freq"] = gw_freq;
                vars["E"] = E;
                vars["B"] = B;
                vars["kg_mass"] = kg_mass;
                vars["spin"] = spin;
                vars["viscosity"] = viscosity;
                vars["temp"] = temp;
                vars["dim"] = dim;
                vars["alpha"] = alpha;
                vars["H0"] = H0;

                try {
                    double value = parser.evaluate(equation, vars);
                    p.ax = value * scale;
                    p.ay = value * scale;
                    p.az = value * scale;

                    p.x += p.vx * sim.dt + 0.5 * p.ax * sim.dt * sim.dt;
                    p.y += p.vy * sim.dt + 0.5 * p.ay * sim.dt * sim.dt;
                    p.z += p.vz * sim.dt + 0.5 * p.az * sim.dt * sim.dt;
                    p.vx += p.ax * sim.dt;
                    p.vy += p.ay * sim.dt;
                    p.vz += p.az * sim.dt;
                } catch (const std::exception& e) {
                    log_message("Failed to evaluate equation '" + equation + "': " + e.what());
                }
            }
        } else {
            #pragma omp parallel
            {
                ExpressionParser local_parser;
                #pragma omp for collapse(3)
                for (int i = 0; i < sim.size; ++i) {
                    for (int j = 0; j < sim.size; ++j) {
                        for (int k = 0; k < sim.size; ++k) {
                            double x = (i - sim.size / 2) * sim.dx;
                            double y = (j - sim.size / 2) * sim.dx;
                            double z = (k - sim.size / 2) * sim.dx;
                            double r = sqrt(x * x + y * y + z * z);

                            std::map<std::string, double> vars = parameters;
                            vars["x"] = x;
                            vars["y"] = y;
                            vars["z"] = z;
                            vars["r"] = r;
                            vars["t"] = sim.t;
                            vars["theta"] = std::acos(z / std::max(r, 1e-6));
                            vars["phi"] = std::atan2(y, x);
                            // Add all member variables to vars
                            vars["scale"] = scale;
                            vars["mass"] = mass;
                            vars["c"] = c;
                            vars["G"] = G;
                            vars["mass1"] = mass1;
                            vars["mass2"] = mass2;
                            vars["wave_speed"] = wave_speed;
                            vars["freq"] = freq;
                            vars["h"] = h;
                            vars["planck_freq"] = planck_freq;
                            vars["a"] = a;
                            vars["T"] = T;
                            vars["M"] = M;
                            vars["radius"] = radius;
                            vars["l"] = l;
                            vars["m"] = m;
                            vars["n"] = n;
                            vars["l_quant"] = l_quant;
                            vars["m_quant"] = m_quant;
                            vars["gw_amp"] = gw_amp;
                            vars["gw_freq"] = gw_freq;
                            vars["E"] = E;
                            vars["B"] = B;
                            vars["kg_mass"] = kg_mass;
                            vars["spin"] = spin;
                            vars["viscosity"] = viscosity;
                            vars["temp"] = temp;
                            vars["dim"] = dim;
                            vars["alpha"] = alpha;
                            vars["H0"] = H0;

                            try {
                                double value = local_parser.evaluate(equation, vars) * scale * vis_scale;
                                sim.ricci[i][j][k].set_diagonal(value);
                                sim.computed_scalar[i][j][k] = value;

                                for (int w = 0; w < sim.size; ++w) {
                                    double w_coord = (w - sim.size / 2) * sim.dx;
                                    sim.scalar_4d[i][j][k][w] = value * (1.0 + 0.3 * sin(w_coord + sim.t));
                                }
                            } catch (const std::exception& e) {
                                log_message("Failed to evaluate equation '" + equation + "': " + e.what());
                            }
                        }
                    }
                }
            }
        }
        sim.t += sim.dt;
    }

    std::string output(const Simulation& sim) const override {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(6);
        oss << "# Equation: " << name() << "\n";
        oss << "# Time: " << sim.t << "\n";
        if (is_particle_based) {
            oss << "# Particles (index,x,y,z,vx,vy,vz,ax,ay,az)\n";
            for (size_t i = 0; i < sim.particles.size(); ++i) {
                const auto& p = sim.particles[i];
                oss << i << "," << p.x << "," << p.y << "," << p.z << ","
                    << p.vx << "," << p.vy << "," << p.vz << ","
                    << p.ax << "," << p.ay << "," << p.az << "\n";
            }
            double total_ke = 0.0;
            for (const auto& p : sim.particles) {
                double speed = std::sqrt(p.vx * p.vx + p.vy * p.vy + p.vz * p.vz);
                total_ke += 0.5 * mass * speed * speed;
            }
            oss << "# Total Kinetic Energy: " << total_ke << "\n";
        } else {
            oss << "# Scalar Field and Gradients (i,j,k,x,y,z,scalar,grad_x,grad_y,grad_z)\n";
            double sum_scalar = 0.0;
            int count = 0;
            for (int i = 0; i < sim.size; ++i) {
                for (int j = 0; j < sim.size; ++j) {
                    for (int k = 0; k < sim.size; ++k) {
                        double x = (i - sim.size / 2) * sim.dx;
                        double y = (j - sim.size / 2) * sim.dx;
                        double z = (k - sim.size / 2) * sim.dx;
                        double scalar = sim.computed_scalar[i][j][k];
                        sum_scalar += scalar;
                        count++;
                        double grad_x = (i > 0 && i < sim.size - 1) ? (sim.computed_scalar[i + 1][j][k] - sim.computed_scalar[i - 1][j][k]) / (2 * sim.dx) : 0.0;
                        double grad_y = (j > 0 && j < sim.size - 1) ? (sim.computed_scalar[i][j + 1][k] - sim.computed_scalar[i][j - 1][k]) / (2 * sim.dx) : 0.0;
                        double grad_z = (k > 0 && k < sim.size - 1) ? (sim.computed_scalar[i][j][k + 1] - sim.computed_scalar[i][j][k - 1]) / (2 * sim.dx) : 0.0;
                        oss << i << "," << j << "," << k << "," << x << "," << y << "," << z << ","
                            << scalar << "," << grad_x << "," << grad_y << "," << grad_z << "\n";
                    }
                }
            }
            oss << "# Average Scalar: " << (count > 0 ? sum_scalar / count : 0.0) << "\n";
            oss << "# Total Scalar (Integral): " << (sum_scalar * sim.dx * sim.dx * sim.dx) << "\n";
        }
        return oss.str();
    }

    std::string name() const override { return eq_name; }

    void setVisScale(double s) { vis_scale = std::max(0.1, std::min(10.0, s)); }
    void setVisColorIntensity(double c) { vis_color_intensity = std::max(0.1, std::min(2.0, c)); }
    double getVisScale() const { return vis_scale; }
    double getVisColorIntensity() const { return vis_color_intensity; }

private:
    std::string eq_name;
    std::string equation;
    ExpressionParser parser;
    std::map<std::string, double> parameters;
    bool is_particle_based;
    bool valid;
    double scale, mass, c, G, mass1, mass2, r, wave_speed, freq, h, planck_freq, a, T, M, radius;
    double l, m, n, l_quant, m_quant, gw_amp, gw_freq, E, B, kg_mass, spin, viscosity, temp, dim, alpha, H0;
    double vis_scale;
    double vis_color_intensity;
};

class TensorEquation : public Equation {
public:
    TensorEquation(const std::string& name, const std::map<std::string, double>& params)
        : eq_name(name), H1D(1.0), D_t_coeff(0.1), eps_5D(0.01), alpha(1.0), I0(0.0), D0(0.1), FiveD0(0.1) {
        for (const auto& param : params) {
            const std::string& key = param.first;
            if (key == "H1D") H1D = param.second;
            else if (key == "D_t") D_t_coeff = param.second;
            else if (key == "eps_5D") eps_5D = param.second;
            else if (key == "alpha") alpha = param.second;
            else if (key == "I0") I0 = param.second;
            else if (key == "D0") D0 = param.second;
            else if (key == "FiveD0") FiveD0 = param.second;
            else {
                log_message("Warning: Unrecognized parameter '" + key + "' in equation '" + name + "'");
            }
        }
    }

    void compute(Simulation& sim) override {
        double D_t = D_t_coeff * sim.t;
        double FiveD_t = FiveD0 * sin(sim.t);

        #pragma omp parallel for collapse(3)
        for (int i = 1; i < sim.size - 1; ++i) {
            for (int j = 1; j < sim.size - 1; ++j) {
                for (int k = 1; k < sim.size - 1; ++k) {
                    Tensor& div = sim.divergence[i][j][k];
                    for (int nu = 0; nu < 3; ++nu) {
                        for (int mu = 0; mu < 3; ++mu) {
                            double dR = 0.0;
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
                    double trace = H1D * alpha * 3.0 + I0 * 3.0 - 4.0 * D_t * 3.0 + eps_5D * FiveD_t * 3.0;
                    rhs.set_diagonal(trace / 3.0);

                    double scalar = 0.0;
                    for (int m = 0; m < 3; ++m) {
                        scalar += div.data[m][m];
                    }
                    sim.computed_scalar[i][j][k] = scalar;

                    for (int m = 0; m < 3; ++m) {
                        for (int n = 0; n < 3; ++n) {
                            sim.ricci[i][j][k].data[m][n] = rhs.data[m][n];
                        }
                    }

                    for (int w = 0; w < sim.size; ++w) {
                        double w_coord = (w - sim.size / 2) * sim.dx;
                        sim.scalar_4d[i][j][k][w] = scalar * (1.0 + 0.1 * sin(w_coord + sim.t));
                    }
                }
            }
        }
        sim.t += sim.dt;
    }

    std::string output(const Simulation& sim) const override {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(6);
        oss << "# Equation: " << name() << "\n";
        oss << "# Time: " << sim.t << "\n";
        oss << "# Scalar Field and Gradients (i,j,k,x,y,z,scalar,grad_x,grad_y,grad_z)\n";
        double sum_scalar = 0.0;
        int count = 0;
        for (int i = 0; i < sim.size; ++i) {
            for (int j = 0; j < sim.size; ++j) {
                for (int k = 0; k < sim.size; ++k) {
                    double x = (i - sim.size / 2) * sim.dx;
                    double y = (j - sim.size / 2) * sim.dx;
                    double z = (k - sim.size / 2) * sim.dx;
                    double scalar = sim.computed_scalar[i][j][k];
                    sum_scalar += scalar;
                    count++;
                    double grad_x = (i > 0 && i < sim.size - 1) ? (sim.computed_scalar[i + 1][j][k] - sim.computed_scalar[i - 1][j][k]) / (2 * sim.dx) : 0.0;
                    double grad_y = (j > 0 && j < sim.size - 1) ? (sim.computed_scalar[i][j + 1][k] - sim.computed_scalar[i][j - 1][k]) / (2 * sim.dx) : 0.0;
                    double grad_z = (k > 0 && k < sim.size - 1) ? (sim.computed_scalar[i][j][k + 1] - sim.computed_scalar[i][j][k - 1]) / (2 * sim.dx) : 0.0;
                    oss << i << "," << j << "," << k << "," << x << "," << y << "," << z << ","
                        << scalar << "," << grad_x << "," << grad_y << "," << grad_z << "\n";
                }
            }
        }
        oss << "# Average Scalar: " << (count > 0 ? sum_scalar / count : 0.0) << "\n";
        oss << "# Total Scalar (Integral): " << (sum_scalar * sim.dx * sim.dx * sim.dx) << "\n";
        return oss.str();
    }

    std::string name() const override { return eq_name; }

private:
    std::string eq_name;
    double H1D, D_t_coeff, eps_5D, alpha, I0, D0, FiveD0;
};

#endif