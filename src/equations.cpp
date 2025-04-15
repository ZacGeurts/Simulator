#include "equations.h"
#include <fstream>
#include <sstream>
#include <ctime>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <vector>
#include <map>
#include <functional>

// Simple expression parser for equations
class ExpressionParser {
public:
    ExpressionParser() {
        functions["sin"] = [](double x) { return std::sin(x); };
        functions["cos"] = [](double x) { return std::cos(x); };
        functions["sqrt"] = [](double x) { return std::sqrt(x); };
        functions["exp"] = [](double x) { return std::exp(x); };
    }

    double evaluate(const std::string& expr, const std::map<std::string, double>& vars) {
        std::string token;
        std::stringstream ss(expr);
        std::vector<std::string> tokens;

        while (ss >> token) {
            tokens.push_back(token);
        }

        std::vector<std::string> stack;
        for (const auto& t : tokens) {
            if (t == "+" || t == "-" || t == "*" || t == "/") {
                stack.push_back(t);
            } else if (functions.count(t)) {
                stack.push_back(t);
            } else {
                if (vars.count(t)) {
                    stack.push_back(std::to_string(vars.at(t)));
                } else {
                    try {
                        std::stod(t);
                        stack.push_back(t);
                    } catch (...) {
                        throw std::runtime_error("Unknown variable or invalid token: " + t);
                    }
                }
            }
        }

        std::vector<double> eval_stack;
        for (const auto& t : stack) {
            if (t == "+") {
                if (eval_stack.size() < 2) throw std::runtime_error("Invalid expression: too few operands for +");
                double b = eval_stack.back(); eval_stack.pop_back();
                double a = eval_stack.back(); eval_stack.pop_back();
                eval_stack.push_back(a + b);
            } else if (t == "-") {
                if (eval_stack.size() < 2) throw std::runtime_error("Invalid expression: too few operands for -");
                double b = eval_stack.back(); eval_stack.pop_back();
                double a = eval_stack.back(); eval_stack.pop_back();
                eval_stack.push_back(a - b);
            } else if (t == "*") {
                if (eval_stack.size() < 2) throw std::runtime_error("Invalid expression: too few operands for *");
                double b = eval_stack.back(); eval_stack.pop_back();
                double a = eval_stack.back(); eval_stack.pop_back();
                eval_stack.push_back(a * b);
            } else if (t == "/") {
                if (eval_stack.size() < 2) throw std::runtime_error("Invalid expression: too few operands for /");
                double b = eval_stack.back(); eval_stack.pop_back();
                double a = eval_stack.back(); eval_stack.pop_back();
                if (b == 0) throw std::runtime_error("Division by zero");
                eval_stack.push_back(a / b);
            } else if (functions.count(t)) {
                if (eval_stack.size() < 1) throw std::runtime_error("Invalid expression: too few operands for " + t);
                double x = eval_stack.back(); eval_stack.pop_back();
                eval_stack.push_back(functions[t](x));
            } else {
                eval_stack.push_back(std::stod(t));
            }
        }

        if (eval_stack.size() != 1) throw std::runtime_error("Invalid expression: incomplete evaluation");
        return eval_stack[0];
    }

private:
    std::map<std::string, std::function<double(double)>> functions;
};

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

// OriginalEquation class
class OriginalEquation : public Equation {
public:
    OriginalEquation(const std::string& name, const std::map<std::string, double>& params)
        : eq_name(name), H1D(1.0), D_t_coeff(0.1), eps_5D(0.01) {
        if (params.count("H1D")) H1D = params.at("H1D");
        if (params.count("D_t")) D_t_coeff = params.at("D_t");
        if (params.count("eps_5D")) eps_5D = params.at("eps_5D");
    }

    void compute(Simulation& sim) override {
        double D_t = D_t_coeff * sin(sim.t);
        double eps = eps_5D * cos(sim.t);

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
                            else
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
    }

    std::string output(const Simulation& sim) const override {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(6);
        oss << "Equation: " << name() << "\n";
        oss << "Time: " << sim.t << "\n";
        oss << "Scalar Field (grid indices i,j,k):\n";
        for (int i = 0; i < sim.size; ++i) {
            for (int j = 0; j < sim.size; ++j) {
                for (int k = 0; k < sim.size; ++k) {
                    oss << i << "," << j << "," << k << ": " << sim.computed_scalar[i][j][k] << "\n";
                }
            }
        }
        return oss.str();
    }

    std::string name() const override { return eq_name; }

private:
    std::string eq_name;
    double H1D, D_t_coeff, eps_5D;
};

// RefinedEquation class
class RefinedEquation : public Equation {
public:
    RefinedEquation(const std::string& name, const std::map<std::string, double>& params)
        : eq_name(name), phi_base(1.0), G(1.0), k4(0.1), k5(0.01) {
        if (params.count("phi")) phi_base = params.at("phi");
        if (params.count("G")) G = params.at("G");
        if (params.count("k4")) k4 = params.at("k4");
        if (params.count("k5")) k5 = params.at("k5");
    }

    void compute(Simulation& sim) override {
        double phi = phi_base + 0.1 * sin(sim.t);

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
    }

    std::string output(const Simulation& sim) const override {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(6);
        oss << "Equation: " << name() << "\n";
        oss << "Time: " << sim.t << "\n";
        oss << "Scalar Field (grid indices i,j,k):\n";
        for (int i = 0; i < sim.size; ++i) {
            for (int j = 0; j < sim.size; ++j) {
                for (int k = 0; k < sim.size; ++k) {
                    oss << i << "," << j << "," << k << ": " << sim.computed_scalar[i][j][k] << "\n";
                }
            }
        }
        return oss.str();
    }

    std::string name() const override { return eq_name; }

private:
    std::string eq_name;
    double phi_base, G, k4, k5;
};

// CustomEquation class
class CustomEquation : public Equation {
public:
    CustomEquation(const std::string& name, const std::vector<std::string>& eqs, const std::map<std::string, double>& params, bool is_particle_based)
        : eq_name(name), equations(eqs), parameters(params), is_particle_based(is_particle_based),
          scale(1.0), mass(1.0), c(299792458.0), G(6.67430e-11), mass1(1.0), mass2(1.0), r(1.0),
          wave_speed(343.0), freq(1.0), h(6.62607015e-34), planck_freq(1.0e15), a(1.0), T(1.0),
          M(1.0), radius(1.0), l(0.0), m(0.0), n(1.0), l_quant(0.0), m_quant(0.0),
          gw_amp(1.0), gw_freq(1.0), E(1.0), B(1.0), kg_mass(1.0), spin(0.5), viscosity(0.01),
          temp(5000.0), dim(10.0), alpha(1.0), H0(70.0) {
        if (params.count("scale")) scale = params.at("scale");
        if (params.count("mass")) mass = params.at("mass");
        if (params.count("c")) c = params.at("c");
        if (params.count("G")) G = params.at("G");
        if (params.count("mass1")) mass1 = params.at("mass1");
        if (params.count("mass2")) mass2 = params.at("mass2");
        if (params.count("r")) r = params.at("r");
        if (params.count("wave_speed")) wave_speed = params.at("wave_speed");
        if (params.count("freq")) freq = params.at("freq");
        if (params.count("h")) h = params.at("h");
        if (params.count("freq")) planck_freq = params.at("freq");
        if (params.count("a")) a = params.at("a");
        if (params.count("T")) T = params.at("T");
        if (params.count("M")) M = params.at("M");
        if (params.count("radius")) radius = params.at("radius");
        if (params.count("l")) l = params.at("l");
        if (params.count("m")) m = params.at("m");
        if (params.count("n")) n = params.at("n");
        if (params.count("l_quant")) l_quant = params.at("l_quant");
        if (params.count("m_quant")) m_quant = params.at("m_quant");
        if (params.count("amp")) gw_amp = params.at("amp");
        if (params.count("freq")) gw_freq = params.at("freq");
        if (params.count("E")) E = params.at("E");
        if (params.count("B")) B = params.at("B");
        if (params.count("mass")) kg_mass = params.at("mass");
        if (params.count("spin")) spin = params.at("spin");
        if (params.count("viscosity")) viscosity = params.at("viscosity");
        if (params.count("T")) temp = params.at("T");
        if (params.count("dim")) dim = params.at("dim");
        if (params.count("alpha")) alpha = params.at("alpha");
        if (params.count("H0")) H0 = params.at("H0");
    }

    void compute(Simulation& sim) override {
        if (is_particle_based) {
            ExpressionParser parser;
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

                for (const auto& eq : equations) {
                    size_t eq_pos = eq.find('=');
                    if (eq_pos == std::string::npos) continue;
                    std::string lhs = eq.substr(0, eq_pos);
                    std::string rhs = eq.substr(eq_pos + 1);
                    while (!lhs.empty() && lhs.back() == ' ') lhs.pop_back();
                    while (!rhs.empty() && rhs[0] == ' ') rhs.erase(0, 1);

                    try {
                        double value = parser.evaluate(rhs, vars);
                        if (lhs == "x") p.x = value;
                        else if (lhs == "y") p.y = value;
                        else if (lhs == "z") p.z = value;
                        else if (lhs == "vx") p.vx = value;
                        else if (lhs == "vy") p.vy = value;
                        else if (lhs == "vz") p.vz = value;
                        else if (lhs == "ax") p.ax = value;
                        else if (lhs == "ay") p.ay = value;
                        else if (lhs == "az") p.az = value;
                        else {
                            parameters[lhs] = value;
                        }
                    } catch (const std::exception& e) {
                        log_message("ERROR: Failed to evaluate equation '" + eq + "': " + e.what());
                    }
                }

                p.x += p.vx * sim.dt + 0.5 * p.ax * sim.dt * sim.dt;
                p.y += p.vy * sim.dt + 0.5 * p.ay * sim.dt * sim.dt;
                p.z += p.vz * sim.dt + 0.5 * p.az * sim.dt * sim.dt;
                p.vx += p.ax * sim.dt;
                p.vy += p.ay * sim.dt;
                p.vz += p.az * sim.dt;
            }
        } else {
            for (int i = 0; i < sim.size; ++i) {
                for (int j = 0; j < sim.size; ++j) {
                    for (int k = 0; k < sim.size; ++k) {
                        Tensor& R = sim.ricci[i][j][k];
                        double value = 0.0;
                        double x = (i - sim.size / 2) * sim.dx;
                        double y = (j - sim.size / 2) * sim.dx;
                        double z = (k - sim.size / 2) * sim.dx;
                        double r = sqrt(x * x + y * y + z * z);

                        if (eq_name == "E_equals_mc2") {
                            value = mass * c * c * (1.0 + 0.1 * sin(sim.t));
                        } else if (eq_name == "Newton_Gravity") {
                            value = (G * mass1 * mass2) / (r * r + 1e-6);
                        } else if (eq_name == "Wave_Equation") {
                            double lambda = wave_speed / freq;
                            value = sin(2 * M_PI * (r / lambda + sim.t));
                        } else if (eq_name == "Planck_Energy") {
                            value = h * planck_freq * (1.0 + 0.1 * sin(sim.t));
                        } else if (eq_name == "Kepler_Third") {
                            value = (T * T) / (a * a * a) * (1.0 + 0.1 * cos(sim.t));
                        } else if (eq_name == "Schwarzschild") {
                            double rs = 2 * G * M / (c * c);
                            value = rs / (r + 1e-6);
                        } else if (eq_name == "Spherical_Harmonic") {
                            double theta = atan2(sqrt(x * x + y * y), z);
                            double phi = atan2(y, x);
                            double Plm = (l == 2 && m == 1) ? sqrt(15.0 / (8 * M_PI)) * sin(theta) * cos(theta) * cos(phi) : 0.0;
                            value = Plm * (1.0 + 0.1 * sin(sim.t));
                        } else if (eq_name == "Hydrogen_Wavefunction") {
                            double r0 = 5.29e-11;
                            double psi = exp(-r / (n * r0)) * (r / r0);
                            value = psi * (1.0 + 0.1 * sin(sim.t));
                        } else if (eq_name == "Gravitational_Wave") {
                            value = gw_amp * sin(2 * M_PI * gw_freq * (sim.t + r));
                        } else if (eq_name == "Electromagnetic_Field") {
                            value = E * E + B * B;
                        } else if (eq_name == "Klein_Gordon") {
                            value = sin(kg_mass * r + sim.t);
                        } else if (eq_name == "Dirac_Field") {
                            value = spin * (1.0 + 0.1 * cos(r + sim.t));
                        } else if (eq_name == "Navier_Stokes") {
                            value = viscosity * (x * x + y * y + z * z);
                        } else if (eq_name == "Blackbody_Radiation") {
                            double lambda = 500e-9;
                            value = (2 * h * c * c) / (pow(lambda, 5) * (exp((h * c) / (lambda * 1.38e-23 * temp)) - 1));
                        } else if (eq_name == "String_Theory") {
                            value = alpha * sin(dim * r + sim.t);
                        } else if (eq_name == "Cosmological_Expansion") {
                            value = H0 * r * (1.0 + 0.1 * sin(sim.t));
                        } else {
                            value = scale * (1.0 + 0.2 * sin(sim.t + r));
                        }

                        R.set_diagonal(value);
                        sim.computed_scalar[i][j][k] = value;

                        for (int w = 0; w < sim.size; ++w) {
                            double w_coord = (w - sim.size / 2) * sim.dx;
                            sim.scalar_4d[i][j][k][w] = value * (1.0 + 0.3 * sin(w_coord + sim.t));
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
        oss << "Equation: " << name() << "\n";
        oss << "Time: " << sim.t << "\n";
        if (is_particle_based) {
            oss << "Particles (index, x, y, z, vx, vy, vz):\n";
            for (size_t i = 0; i < sim.particles.size(); ++i) {
                const auto& p = sim.particles[i];
                oss << i << ": " << p.x << ", " << p.y << ", " << p.z << ", "
                    << p.vx << ", " << p.vy << ", " << p.vz << "\n";
            }
        } else {
            oss << "Scalar Field (grid indices i,j,k):\n";
            for (int i = 0; i < sim.size; ++i) {
                for (int j = 0; j < sim.size; ++j) {
                    for (int k = 0; k < sim.size; ++k) {
                        oss << i << "," << j << "," << k << ": " << sim.computed_scalar[i][j][k] << "\n";
                    }
                }
            }
        }
        return oss.str();
    }

    std::string name() const override { return eq_name; }

private:
    std::string eq_name;
    std::vector<std::string> equations;
    mutable std::map<std::string, double> parameters;
    bool is_particle_based;
    double scale, mass, c, G, mass1, mass2, r, wave_speed, freq, h, planck_freq, a, T, M, radius;
    double l, m, n, l_quant, m_quant, gw_amp, gw_freq, E, B, kg_mass, spin, viscosity, temp, dim, alpha, H0;
};

// Simulation class implementation
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
        std::cerr << "Failed to open " << filename << std::endl;
        log_message("ERROR: Failed to open " + filename);
        return;
    }
    log_message("Opened " + filename);

    std::string line;
    int eq_count = 1;
    std::vector<std::string> current_eq_lines;
    std::string current_eq_name;
    std::map<std::string, double> params;
    bool in_equation_block = false;

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') {
            log_message("Skipped line: " + line);
            continue;
        }
        log_message("Processing line: " + line);

        if (line.find("param") == 0) {
            std::istringstream iss(line);
            std::string param_keyword, token;
            iss >> param_keyword;
            while (iss >> token) {
                size_t eq_pos = token.find('=');
                if (eq_pos != std::string::npos) {
                    std::string key = token.substr(0, eq_pos);
                    std::string value_str = token.substr(eq_pos + 1);
                    try {
                        double value = std::stod(value_str);
                        if (key == "dt") dt = value;
                        params[key] = value;
                        log_message("Parsed parameter: " + key + "=" + std::to_string(value));
                    } catch (const std::exception& e) {
                        log_message("ERROR: Invalid parameter value: " + token + " (" + e.what() + ")");
                    }
                }
            }
            continue;
        }

        if (line.find("equation") == 0) {
            if (in_equation_block && !current_eq_name.empty() && !current_eq_lines.empty()) {
                equations.push_back(new CustomEquation(current_eq_name, current_eq_lines, params, true));
                log_message("Added particle-based equation: " + current_eq_name);
                compute();
                save_output();
                switch_equation(equations.size() - 1);
            }
            current_eq_lines.clear();
            params.clear();
            std::istringstream iss(line);
            std::string keyword;
            iss >> keyword >> current_eq_name;
            in_equation_block = true;
            continue;
        }

        if (in_equation_block) {
            size_t eq_pos = line.find('=');
            if (eq_pos != std::string::npos && line.find("equation") == std::string::npos) {
                current_eq_lines.push_back(line);
            } else {
                std::istringstream iss(line);
                std::string token;
                while (iss >> token) {
                    eq_pos = token.find('=');
                    if (eq_pos != std::string::npos) {
                        std::string key = token.substr(0, eq_pos);
                        std::string value_str = token.substr(eq_pos + 1);
                        try {
                            double value = std::stod(value_str);
                            params[key] = value;
                            log_message("Parsed parameter for " + current_eq_name + ": " + key + "=" + std::to_string(value));
                        } catch (const std::exception& e) {
                            log_message("ERROR: Invalid parameter value: " + token + " (" + e.what() + ")");
                        }
                    }
                }
            }
            continue;
        }

        size_t param_start = std::string::npos;
        for (size_t i = 0; i < line.size(); ++i) {
            if (line[i] == '=' && (i == 0 || line[i-1] != ' ') && (i+1 < line.size() && std::isdigit(line[i+1]))) {
                param_start = i;
                break;
            }
        }
        std::string equation_str;
        std::string param_str;
        if (param_start != std::string::npos) {
            equation_str = line.substr(0, param_start);
            while (equation_str.back() == ' ') equation_str.pop_back();
            param_str = line.substr(param_start);
        } else {
            equation_str = line;
        }
        log_message("Equation string: " + equation_str);

        std::map<std::string, double> local_params;
        if (!param_str.empty()) {
            std::istringstream iss(param_str);
            std::string token;
            while (iss >> token) {
                size_t eq_pos = token.find('=');
                if (eq_pos != std::string::npos) {
                    std::string key = token.substr(0, eq_pos);
                    std::string value_str = token.substr(eq_pos + 1);
                    try {
                        if (!value_str.empty()) {
                            double value = std::stod(value_str);
                            local_params[key] = value;
                            log_message("Parsed parameter: " + key + "=" + std::to_string(value));
                        }
                    } catch (const std::exception& e) {
                        log_message("ERROR: Invalid parameter value: " + token + " (" + e.what() + ")");
                    }
                }
            }
        }

        std::string name = "Eq_" + std::to_string(eq_count++);
        if (equation_str.find("∇") != std::string::npos || equation_str.find("\\nabla") != std::string::npos) {
            equations.push_back(new OriginalEquation(equation_str, local_params));
            log_message("Added original equation: " + equation_str);
            compute();
            save_output();
            switch_equation(equations.size() - 1);
        } else {
            std::istringstream iss2(equation_str);
            std::string type;
            iss2 >> name >> type;
            if (type == "refined") {
                equations.push_back(new RefinedEquation(name, local_params));
                log_message("Added refined equation: " + name);
                compute();
                save_output();
                switch_equation(equations.size() - 1);
            } else if (type == "custom") {
                equations.push_back(new CustomEquation(name, {}, local_params, false));
                log_message("Added custom equation: " + name);
                compute();
                save_output();
                switch_equation(equations.size() - 1);
            } else {
                log_message("WARNING: Unknown equation type, skipping: " + type);
            }
        }
    }

    if (in_equation_block && !current_eq_name.empty() && !current_eq_lines.empty()) {
        equations.push_back(new CustomEquation(current_eq_name, current_eq_lines, params, true));
        log_message("Added particle-based equation: " + current_eq_name);
        compute();
        save_output();
        switch_equation(equations.size() - 1);
    }

    file.close();
    if (equations.empty()) {
        std::cerr << "No equations loaded" << std::endl;
        log_message("ERROR: No equations loaded");
    } else {
        switch_equation(0);
        log_message("Loaded " + std::to_string(equations.size()) + " equations");
    }
}

void Simulation::compute() {
    if (!equations.empty()) {
        equations[current_equation]->compute(*this);
    }
}

void Simulation::save_output() {
    if (equations.empty()) return;
    std::string eq_name = equations[current_equation]->name();
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
        file << equations[current_equation]->output(*this);
        file.close();
        output_written.insert(eq_name);
        log_message("Wrote output: output/equation_" + safe_name + ".txt");
    } else {
        log_message("ERROR: Failed to write output: output/equation_" + safe_name + ".txt");
    }
}

void Simulation::switch_equation(int index) {
    current_equation = index;
    t = 0.0;
    compute();
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
    if (equations.empty()) return;
    current_equation = (current_equation - 1 + equations.size()) % equations.size();
    t = 0.0;
    compute();
    save_output();
}

void Simulation::next_equation() {
    if (equations.empty()) return;
    current_equation = (current_equation + 1) % equations.size();
    t = 0.0;
    compute();
    save_output();
}