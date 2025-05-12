#include "equations.h"
#include <fstream>
#include <sstream>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <mutex>

// Thread-safe buffered logger
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
        if (log_file.is_open()) {
            for (const auto& msg : buffer) {
                log_file << msg << std::endl;
                std::cout << msg << std::endl;
            }
            buffer.clear();
            log_file.close();
        }
    }
};

static Logger logger;

// Tensor implementation
Tensor::Tensor() {
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            data[i][j] = 0.0;
}

void Tensor::set_diagonal(double value) {
    data[0][0] = data[1][1] = data[2][2] = value;
}

// Particle implementation
Particle::Particle() : x(0.0), y(0.0), z(0.0), vx(0.0), vy(0.0), vz(0.0), ax(0.0), ay(0.0), az(0.0) {}

// ExpressionParser implementation
ExpressionParser::ExpressionParser() {
    functions["sin"] = [](double x) { return std::sin(x); };
    functions["cos"] = [](double x) { return std::cos(x); };
    functions["sqrt"] = [](double x) { return std::sqrt(x); };
    functions["exp"] = [](double x) { return std::exp(x); };
    functions["log"] = [](double x) { return std::log(x); };
    functions["tan"] = [](double x) { return std::tan(x); };
    functions["acos"] = [](double x) { return std::acos(x); };
    functions["asin"] = [](double x) { return std::asin(x); };
    functions["atan"] = [](double x) { return std::atan(x); };
    binary_functions["max"] = [](double x, double y) { return std::max(x, y); };
    binary_functions["min"] = [](double x, double y) { return std::min(x, y); };
    binary_functions["atan2"] = [](double y, double x) { return std::atan2(y, x); };
}

double ExpressionParser::evaluate(const std::string& expr, const std::map<std::string, double>& vars) {
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

bool ExpressionParser::validate(const std::string& expr) const {
    for (char c : expr) {
        if (!std::isalnum(c) && c != '_' && c != '+' && c != '-' && c != '*' && c != '/' &&
            c != '^' && c != '(' && c != ')' && c != '.' && c != ' ' && c != '=' && c != ',' &&
            c != 'e') {
            return false;
        }
    }
    return true;
}

std::string ExpressionParser::preprocess(const std::string& expr) {
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
        last_was_digit = std::isdigit(c) || c == '.' || c == 'e';
        last_was_var = std::isalpha(c);
    }
    return result;
}

void ExpressionParser::skipWhitespace() {
    while (pos < expression.length() && std::isspace(expression[pos])) {
        ++pos;
    }
}

double ExpressionParser::parseExpression() {
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

double ExpressionParser::parseTerm() {
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

double ExpressionParser::parseFactor() {
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

double ExpressionParser::parseNumber() {
    skipWhitespace();
    size_t start = pos;
    while (pos < expression.length() && (std::isdigit(expression[pos]) || expression[pos] == '.' || expression[pos] == 'e' || expression[pos] == '+' || expression[pos] == '-')) {
        ++pos;
    }
    std::string num_str = expression.substr(start, pos - start);
    try {
        size_t idx;
        double value = std::stod(num_str, &idx);
        if (idx != num_str.length()) {
            throw std::runtime_error("Invalid number format: " + num_str);
        }
        return value;
    } catch (...) {
        throw std::runtime_error("Invalid number: " + num_str);
    }
}

std::string ExpressionParser::parseToken() {
    skipWhitespace();
    size_t start = pos;
    while (pos < expression.length() && (std::isalnum(expression[pos]) || expression[pos] == '_')) {
        ++pos;
    }
    if (start == pos) throw std::runtime_error("Expected token at position " + std::to_string(pos));
    return expression.substr(start, pos - start);
}

// Simulation implementation
Simulation::Simulation() : t(0.0), dt(0.016), current_equation(0), size(0), dx(0.0) {}

void Simulation::initialize(int num_particles, int grid_size, double grid_dx, double dt_val) {
    dt = dt_val;
    t = 0.0;
    size = grid_size;
    dx = grid_dx;

    particles.resize(num_particles);
    for (int i = 0; i < num_particles; ++i) {
        particles[i].x = (i - num_particles / 2) * 1.0;
        particles[i].y = 0.0;
        particles[i].z = 0.0;
        particles[i].vx = 0.0;
        particles[i].vy = 0.0;
        particles[i].vz = 0.0;
    }

    ricci.resize(size, std::vector<std::vector<Tensor>>(size, std::vector<Tensor>(size)));
    divergence.resize(size, std::vector<std::vector<Tensor>>(size, std::vector<Tensor>(size)));
    scalar.resize(size, std::vector<std::vector<double>>(size, std::vector<double>(size, 0.0)));
    computed_scalar.resize(size, std::vector<std::vector<double>>(size, std::vector<double>(size, 0.0)));
    scalar_4d.resize(size, std::vector<std::vector<std::vector<double>>>(size, std::vector<std::vector<double>>(size, std::vector<double>(size, 0.0))));

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            for (int k = 0; k < size; ++k) {
                double x = (i - size / 2) * dx;
                double y = (j - size / 2) * dx;
                double z = (k - size / 2) * dx;
                double r2 = x * x + y * y + z * z;
                ricci[i][j][k].set_diagonal(exp(-r2));
                scalar[i][j][k] = exp(-r2);
                computed_scalar[i][j][k] = scalar[i][j][k];
                for (int w = 0; w < size; ++w) {
                    double w_coord = (w - size / 2) * dx;
                    scalar_4d[i][j][k][w] = exp(-r2) * (1.0 + 0.1 * sin(w_coord));
                }
            }
        }
    }
    log_message("Initialized simulation: num_particles=" + std::to_string(num_particles) +
                ", grid_size=" + std::to_string(size) + ", dx=" + std::to_string(dx) + ", dt=" + std::to_string(dt));
}

void Simulation::load_equations(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        log_message("ERROR: Failed to open " + filename);
        std::cerr << "Failed to open " << filename << std::endl;
        return;
    }
    log_message("Opened " + filename);

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        log_message("Processing line: " + line);

        std::istringstream iss(line);
        std::string name, type, equation;
        if (!(iss >> name >> type)) {
            log_message("ERROR: Invalid line format, expected 'name type': " + line);
            continue;
        }
        std::getline(iss, equation);
        equation.erase(0, equation.find_first_not_of(" \t"));

        std::map<std::string, double> params;
        std::map<std::string, std::string> derived_params;
        size_t param_start = equation.find_first_of(" \t", equation.find_first_not_of(" \t"));
        if (param_start != std::string::npos) {
            std::string param_str = equation.substr(param_start);
            equation = equation.substr(0, param_start);
            equation.erase(equation.find_last_not_of(" \t") + 1);
            std::istringstream param_iss(param_str);
            std::string token;
            while (param_iss >> token) {
                size_t eq_pos = token.find('=');
                if (eq_pos != std::string::npos) {
                    std::string key = token.substr(0, eq_pos);
                    std::string value_str = token.substr(eq_pos + 1);
                    try {
                        // Check if the value is a derived expression
                        if (value_str.find('(') != std::string::npos || value_str.find('/') != std::string::npos) {
                            derived_params[key] = value_str;
                            params[key] = 0.0; // Placeholder
                            log_message("Recognized derived parameter: " + key + "=" + value_str);
                            continue;
                        }
                        // Parse numeric value
                        size_t idx;
                        double value = std::stod(value_str, &idx);
                        if (idx != value_str.length()) {
                            throw std::invalid_argument("Invalid numeric format");
                        }
                        if (key == "dt") dt = value;
                        params[key] = value;
                        log_message("Parsed parameter: " + key + "=" + std::to_string(value));
                    } catch (const std::exception& e) {
                        log_message("ERROR: Invalid parameter value: " + token + ", " + e.what());
                    }
                }
            }
        }

        if (type == "custom" || type == "particle") {
            size_t eq_pos = equation.find('=');
            if (eq_pos != std::string::npos) {
                equation = equation.substr(eq_pos + 1);
                equation.erase(0, equation.find_first_not_of(" \t"));
                equation.erase(equation.find_last_not_of(" \t") + 1);
            }
        }

        try {
            Equation* eq = nullptr;
            if (type == "custom") {
                eq = new CustomEquation(name, equation, params, derived_params, false);
            } else if (type == "particle") {
                eq = new CustomEquation(name, equation, params, derived_params, true);
            } else if (type == "tensor") {
                eq = new TensorEquation(name, equation, params, derived_params);
            } else if (type == "quantum") {
                eq = new QuantumEquation(name, params);
            } else if (type == "original") {
                eq = new OriginalEquation(name, params);
            } else if (type == "refined") {
                eq = new RefinedEquation(name, params);
            } else {
                log_message("WARNING: Unknown equation type: " + type);
                continue;
            }

            if (eq->is_valid()) {
                equations.push_back(eq);
                log_message("Added equation: " + name + ", type: " + type);
            } else {
                log_message("WARNING: Skipping invalid equation: " + name);
                delete eq;
            }
        } catch (const std::exception& e) {
            log_message("ERROR: Failed to create equation '" + name + "': " + e.what());
        }
    }

    file.close();
    logger.flush();
    if (equations.empty()) {
        log_message("ERROR: No valid equations loaded");
        std::cerr << "No valid equations loaded" << std::endl;
    } else {
        switch_equation(0);
        log_message("Loaded " + std::to_string(equations.size()) + " equations");
    }
}

void Simulation::compute() {
    if (equations.empty()) {
        log_message("ERROR: No equations to compute");
        return;
    }
    if (current_equation < 0 || current_equation >= static_cast<int>(equations.size())) {
        log_message("ERROR: Invalid current_equation index: " + std::to_string(current_equation));
        return;
    }
    log_message("Computing equation: " + equations[current_equation]->name());
    try {
        equations[current_equation]->compute(*this);
    } catch (const std::exception& e) {
        log_message("ERROR: Compute failed for equation '" + equations[current_equation]->name() + "': " + e.what());
    }
    logger.flush();
}

void Simulation::save_output() {
    if (equations.empty()) {
        log_message("ERROR: No equations to save output");
        return;
    }
    if (current_equation < 0 || current_equation >= static_cast<int>(equations.size())) {
        log_message("ERROR: Invalid current_equation index: " + std::to_string(current_equation));
        return;
    }
    std::string eq_name = equations[current_equation]->name();
    log_message("Saving output for equation: " + eq_name);
    if (output_written.count(eq_name)) return;

    std::string safe_name = eq_name;
    std::string temp_name;
    bool last_was_underscore = false;
    for (char c : safe_name) {
        if (std::isalnum(c) || c == '(' || c == ')' || c == '=') {
            temp_name += c;
            last_was_underscore = false;
        } else if (c == ' ') {
            if (!last_was_underscore) {
                temp_name += '_';
                last_was_underscore = true;
            }
        } else {
            if (!last_was_underscore) {
                temp_name += '_';
                last_was_underscore = true;
            }
        }
    }
    safe_name = temp_name;
    if (!safe_name.empty() && safe_name.back() == '_') {
        safe_name.pop_back();
    }

    std::ofstream file("output/equation_" + safe_name + ".txt");
    if (file.is_open()) {
        file << "# Equation output at time t=" << t << " (units: meters, seconds)\n";
        try {
            file << equations[current_equation]->output(*this);
        } catch (const std::exception& e) {
            log_message("ERROR: Output failed for equation '" + eq_name + "': " + e.what());
        }
        file.close();
        output_written.insert(eq_name);
        log_message("Wrote output: output/equation_" + safe_name + ".txt");
    } else {
        log_message("ERROR: Failed to write output: output/equation_" + safe_name + ".txt");
    }
    logger.flush();
}

void Simulation::switch_equation(int index) {
    if (index >= 0 && index < static_cast<int>(equations.size())) {
        current_equation = index;
        t = 0.0;
        log_message("Switched to equation: " + equations[current_equation]->name());
        compute();
        save_output();
    } else {
        log_message("ERROR: Invalid equation index: " + std::to_string(index));
    }
}

void Simulation::setup() {
    initialize(100, 10, 0.1, 0.016);
    load_equations("equations.txt");
    compute();
    save_output();
}

void Simulation::step() {
    compute();
    save_output();
}

void Simulation::prev_equation() {
    if (equations.empty()) {
        log_message("ERROR: No equations to switch");
        return;
    }
    current_equation = (current_equation - 1 + equations.size()) % equations.size();
    t = 0.0;
    log_message("Switched to previous equation: " + equations[current_equation]->name());
    compute();
    save_output();
}

void Simulation::next_equation() {
    if (equations.empty()) {
        log_message("ERROR: No equations to switch");
        return;
    }
    current_equation = (current_equation + 1) % equations.size();
    t = 0.0;
    log_message("Switched to next equation: " + equations[current_equation]->name());
    compute();
    save_output();
}

// CustomEquation implementation
CustomEquation::CustomEquation(const std::string& name, const std::string& eq, const std::map<std::string, double>& params,
                             const std::map<std::string, std::string>& derived_params, bool is_particle_based)
    : eq_name(name), equation(eq), parameters(params), derived_parameters(derived_params),
      is_particle_based(is_particle_based), valid(parser.validate(eq)),
      scale(1.0), mass(1.0), c(299792458.0), G(6.67430e-11), mass1(1.0), mass2(1.0), r(1.0),
      wave_speed(343.0), freq(1.0), h(6.62607015e-34), planck_freq(1.0e15), a(1.0), T(1.0),
      M(1.0), radius(1.0), l(0.0), m(0.0), n(1.0), l_quant(0.0), m_quant(0.0),
      gw_amp(1.0), gw_freq(1.0), E(1.0), B(1.0), kg_mass(1.0), spin(0.5), viscosity(0.01),
      temp(5000.0), dim(10.0), alpha(1.0), H0(70.0), H1D(1.0), D_t(0.1), eps_5D(0.01),
      I0(0.0), D0(0.1), FiveD0(0.1), k(1.0), w(1.0), a0(5.29e-11), vis_scale(1.0), vis_color_intensity(1.0) {
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
            log_message("INFO: Accepted parameter: " + key + "=" + std::to_string(param.second));
        }
    }
    if (!valid) {
        log_message("Invalid equation syntax: " + eq);
    }
}

bool CustomEquation::is_valid() const {
    return valid;
}

void CustomEquation::compute(Simulation& sim) {
    if (!valid) {
        log_message("ERROR: Cannot compute invalid equation: " + eq_name);
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
            vars["theta"] = vars["r"] > 1e-6 ? std::acos(p.z / vars["r"]) : 0.0;
            vars["phi"] = std::atan2(p.y, p.x);
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
            vars["H1D"] = H1D;
            vars["D_t"] = D_t;
            vars["eps_5D"] = eps_5D;
            vars["I0"] = I0;
            vars["D0"] = D0;
            vars["FiveD0"] = FiveD0;
            vars["k"] = k;
            vars["w"] = w;
            vars["a0"] = a0;

            for (const auto& dp : derived_parameters) {
                try {
                    vars[dp.first] = parser.evaluate(dp.second, vars);
                } catch (const std::exception& e) {
                    log_message("ERROR: Failed to evaluate derived parameter '" + dp.first + "': " + e.what());
                    vars[dp.first] = 0.0;
                }
            }

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
                log_message("ERROR: Failed to evaluate equation '" + equation + "': " + e.what());
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
                        vars["theta"] = r > 1e-6 ? std::acos(z / r) : 0.0;
                        vars["phi"] = std::atan2(y, x);
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
                        vars["H1D"] = H1D;
                        vars["D_t"] = D_t;
                        vars["eps_5D"] = eps_5D;
                        vars["I0"] = I0;
                        vars["D0"] = D0;
                        vars["FiveD0"] = FiveD0;
                        vars["k"] = k;
                        vars["w"] = w;
                        vars["a0"] = a0;

                        for (const auto& dp : derived_parameters) {
                            try {
                                vars[dp.first] = local_parser.evaluate(dp.second, vars);
                            } catch (const std::exception& e) {
                                log_message("ERROR: Failed to evaluate derived parameter '" + dp.first + "': " + e.what());
                                vars[dp.first] = 0.0;
                            }
                        }

                        try {
                            double value = local_parser.evaluate(equation, vars) * scale * vis_scale;
                            sim.ricci[i][j][k].set_diagonal(value);
                            sim.computed_scalar[i][j][k] = value;

                            for (int w = 0; w < sim.size; ++w) {
                                double w_coord = (w - sim.size / 2) * sim.dx;
                                sim.scalar_4d[i][j][k][w] = value * (1.0 + 0.3 * sin(w_coord + sim.t));
                            }
                        } catch (const std::exception& e) {
                            log_message("ERROR: Failed to evaluate equation '" + equation + "': " + e.what());
                        }
                    }
                }
            }
        }
    }
    sim.t += sim.dt;
    logger.flush();
}

std::string CustomEquation::output(const Simulation& sim) const {
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
        double total_Ke = 0.0;
        for (const auto& p : sim.particles) {
            double speed = std::sqrt(p.vx * p.vx + p.vy * p.vy + p.vz * p.vz);
            total_Ke += 0.5 * mass * speed * speed;
        }
        oss << "# Total Kinetic Energy: " << total_Ke << "\n";
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

std::string CustomEquation::name() const {
    return eq_name;
}

void CustomEquation::setVisScale(double s) {
    vis_scale = std::max(0.1, std::min(10.0, s));
}

double CustomEquation::getVisScale() const {
    return vis_scale;
}

void CustomEquation::setVisColorIntensity(double c) {
    vis_color_intensity = std::max(0.1, std::min(2.0, c));
}

double CustomEquation::getVisColorIntensity() const {
    return vis_color_intensity;
}

// TensorEquation implementation
TensorEquation::TensorEquation(const std::string& name, const std::string& eq, const std::map<std::string, double>& params,
                             const std::map<std::string, std::string>& derived_params)
    : eq_name(name), equation(eq), parameters(params), derived_parameters(derived_params),
      H1D(1.0), D_t_coeff(0.1), eps_5D(0.01), alpha(1.0), I0(0.0), D0(0.1), FiveD0(0.1) {
    for (const auto& param : params) {
        const std::string& key = param.first;
        if (key == "H1D") H1D = param.second;
        else if (key == "D_t") D_t_coeff = param.second;
        else if (key == "eps_5D") eps_5D = param.second;
        else if (key == "alpha") alpha = param.second;
        else if (key == "I0") I0 = param.second;
        else if (key == "D0") D0 = param.second;
        else if (key == "FiveD0") FiveD0 = param.second;
        else if (key == "vis_scale" || key == "vis_color_intensity") {
            continue;
        } else {
            log_message("INFO: Accepted parameter: " + key + "=" + std::to_string(param.second));
        }
    }
}

void TensorEquation::compute(Simulation& sim) {
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
                    vars["theta"] = r > 1e-6 ? std::acos(z / r) : 0.0;
                    vars["phi"] = std::atan2(y, x);
                    vars["H1D"] = H1D;
                    vars["D_t"] = D_t_coeff;
                    vars["eps_5D"] = eps_5D;
                    vars["alpha"] = alpha;
                    vars["I0"] = I0;
                    vars["D0"] = D0;
                    vars["FiveD0"] = FiveD0;

                    for (const auto& dp : derived_parameters) {
                        try {
                            vars[dp.first] = local_parser.evaluate(dp.second, vars);
                        } catch (const std::exception& e) {
                            log_message("ERROR: Failed to evaluate derived parameter '" + dp.first + "': " + e.what());
                            vars[dp.first] = 0.0;
                        }
                    }

                    try {
                        double value = local_parser.evaluate(equation, vars);
                        sim.ricci[i][j][k].set_diagonal(value);
                        sim.computed_scalar[i][j][k] = value;

                        // Compute divergence for tensor field
                        Tensor& div = sim.divergence[i][j][k];
                        if (i > 0 && i < sim.size - 1 && j > 0 && j < sim.size - 1 && k > 0 && k < sim.size - 1) {
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
                        }

                        for (int w = 0; w < sim.size; ++w) {
                            double w_coord = (w - sim.size / 2) * sim.dx;
                            sim.scalar_4d[i][j][k][w] = value * (1.0 + 0.1 * sin(w_coord + sim.t));
                        }
                    } catch (const std::exception& e) {
                        log_message("ERROR: Failed to evaluate equation '" + equation + "': " + e.what());
                    }
                }
            }
        }
    }
    sim.t += sim.dt;
    logger.flush();
}

std::string TensorEquation::output(const Simulation& sim) const {
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

std::string TensorEquation::name() const {
    return eq_name;
}

// QuantumEquation implementation
QuantumEquation::QuantumEquation(const std::string& name, const std::map<std::string, double>& params)
    : eq_name(name), hbar(1.05457182e-34), m(9.1093837e-31), a0(5.29e-11) {
    for (const auto& param : params) {
        const std::string& key = param.first;
        if (key == "hbar") hbar = param.second;
        else if (key == "m") m = param.second;
        else if (key == "a0") a0 = param.second;
        else if (key == "vis_scale" || key == "vis_color_intensity") {
            continue;
        } else {
            log_message("INFO: Accepted parameter: " + key + "=" + std::to_string(param.second));
        }
    }
}

void QuantumEquation::compute(Simulation& sim) {
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < sim.size; ++i) {
        for (int j = 0; j < sim.size; ++j) {
            for (int k = 0; k < sim.size; ++k) {
                double x = (i - sim.size / 2) * sim.dx;
                double y = (j - sim.size / 2) * sim.dx;
                double z = (k - sim.size / 2) * sim.dx;
                double r = sqrt(x * x + y * y + z * z) + 1e-6;
                double psi = exp(-r / a0) / (sqrt(M_PI) * a0 * a0);
                sim.computed_scalar[i][j][k] = psi * psi;
                sim.ricci[i][j][k].set_diagonal(psi * psi);

                for (int w = 0; w < sim.size; ++w) {
                    double w_coord = (w - sim.size / 2) * sim.dx;
                    sim.scalar_4d[i][j][k][w] = psi * psi * (1.0 + 0.1 * cos(w_coord + sim.t));
                }
            }
        }
    }
    sim.t += sim.dt;
    logger.flush();
}

std::string QuantumEquation::output(const Simulation& sim) const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6);
    oss << "# Equation: " << name() << "\n";
    oss << "# Time: " << sim.t << "\n";
    oss << "# Scalar Field (Probability Density) and Gradients (i,j,k,x,y,z,scalar,grad_x,grad_y,grad_z)\n";
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

std::string QuantumEquation::name() const {
    return eq_name;
}

// OriginalEquation implementation
OriginalEquation::OriginalEquation(const std::string& name, const std::map<std::string, double>& params)
    : eq_name(name), H1D(1.0), D_t_coeff(0.1), eps_5D(0.01) {
    for (const auto& param : params) {
        const std::string& key = param.first;
        if (key == "H1D") H1D = param.second;
        else if (key == "D_t") D_t_coeff = param.second;
        else if (key == "eps_5D") eps_5D = param.second;
        else if (key == "vis_scale" || key == "vis_color_intensity") {
            continue;
        } else {
            log_message("INFO: Accepted parameter: " + key + "=" + std::to_string(param.second));
        }
    }
}

void OriginalEquation::compute(Simulation& sim) {
    double D_t = D_t_coeff * sin(sim.t);
    double eps = eps_5D * cos(sim.t);

    #pragma omp parallel for collapse(3)
    for (int i = 1; i < sim.size - 1; ++i) {
        for (int j = 1; j < sim.size - 1; ++j) {
            for (int k = 1; k < sim.size - 1; ++k) {
                Tensor& div = sim.divergence[i][j][k];
                for (int nu = 0; nu < 3; ++nu) {
                    for (int mu = 0; mu < 3; ++mu) {
                        double dR = 0.0;
                        if (mu == 0)
                            dR = (sim.ricci[i + 1][j][k].data[mu][nu] - sim.ricci[i - 1][j][k].data[mu][nu]) / (2 * sim.dx);
                        else if (mu == 1)
                            dR = (sim.ricci[i][j + 1][k].data[mu][nu] - sim.ricci[i][j - 1][k].data[mu][nu]) / (2 * sim.dx);
                        else if (mu == 2)
                            dR = (sim.ricci[i][j][k + 1].data[mu][nu] - sim.ricci[i][j][k - 1].data[mu][nu]) / (2 * sim.dx);
                        div.data[mu][nu] = dR;
                    }
                }

                Tensor rhs;
                rhs.set_diagonal(H1D);
                for (int m = 0; m < 3; ++m)
                    rhs.data[m][m] += sim.ricci[i][j][k].data[m][m];
                for (int m = 0; m < 3; ++m)
                    rhs.data[m][m] -= 4 * D_t;
                for (int m = 0; m < 3; ++m)
                    rhs.data[m][m] += eps;

                sim.computed_scalar[i][j][k] = 0.0;
                for (int m = 0; m < 3; ++m)
                    sim.computed_scalar[i][j][k] += div.data[m][m] * div.data[m][m];
                sim.computed_scalar[i][j][k] = sqrt(sim.computed_scalar[i][j][k]);

                for (int w = 0; w < sim.size; ++w) {
                    double w_coord = (w - sim.size / 2) * sim.dx;
                    sim.scalar_4d[i][j][k][w] = sim.computed_scalar[i][j][k] * (1.0 + 0.1 * sin(w_coord + sim.t));
                }
            }
        }
    }
    sim.t += sim.dt;
    logger.flush();
}

std::string OriginalEquation::output(const Simulation& sim) const {
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

std::string OriginalEquation::name() const {
    return eq_name;
}

// RefinedEquation implementation
RefinedEquation::RefinedEquation(const std::string& name, const std::map<std::string, double>& params)
    : eq_name(name), phi_base(1.0), G(1.0), k4(0.1), k5(0.01) {
    for (const auto& param : params) {
        const std::string& key = param.first;
        if (key == "phi") phi_base = param.second;
        else if (key == "G") G = param.second;
        else if (key == "k4") k4 = param.second;
        else if (key == "k5") k5 = param.second;
        else if (key == "vis_scale" || key == "vis_color_intensity") {
            continue;
        } else {
            log_message("INFO: Accepted parameter: " + key + "=" + std::to_string(param.second));
        }
    }
}

void RefinedEquation::compute(Simulation& sim) {
    double phi = phi_base + 0.1 * sin(sim.t);

    #pragma omp parallel for collapse(3)
    for (int i = 0; i < sim.size; ++i) {
        for (int j = 0; j < sim.size; ++j) {
            for (int k = 0; k < sim.size; ++k) {
                Tensor& R = sim.ricci[i][j][k];
                Tensor rhs;
                rhs.set_diagonal(phi);
                for (int m = 0; m < 3; ++m) {
                    double x = (i - sim.size / 2) * sim.dx;
                    double y = (j - sim.size / 2) * sim.dx;
                    double z = (k - sim.size / 2) * sim.dx;
                    double r2 = x * x + y * y + z * z;
                    rhs.data[m][m] += 8 * M_PI * G * exp(-r2);
                }
                for (int m = 0; m < 3; ++m)
                    rhs.data[m][m] += k4 * 0.1 * sin(sim.t);
                for (int m = 0; m < 3; ++m)
                    rhs.data[m][m] += k5 * 0.01 * cos(sim.t);

                for (int m = 0; m < 3; ++m)
                    for (int n = 0; n < 3; ++n)
                        R.data[m][n] = rhs.data[m][n];

                sim.computed_scalar[i][j][k] = 0.0;
                for (int m = 0; m < 3; ++m)
                    sim.computed_scalar[i][j][k] += R.data[m][m];

                for (int w = 0; w < sim.size; ++w) {
                    double w_coord = (w - sim.size / 2) * sim.dx;
                    sim.scalar_4d[i][j][k][w] = sim.computed_scalar[i][j][k] * (1.0 + 0.2 * cos(w_coord + sim.t));
                }
            }
        }
    }
    sim.t += sim.dt;
    logger.flush();
}

std::string RefinedEquation::output(const Simulation& sim) const {
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

std::string RefinedEquation::name() const {
    return eq_name;
}