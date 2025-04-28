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
#include <cctype>

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
    double phi_base, G, k4, k5;
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
        log_message("ERROR: Failed to open " + filename);
        std::cerr << "Failed to open " << filename << std::endl;
        return;
    }
    log_message("Opened " + filename);

    std::string line;
    int eq_count = 1;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') {
            //log_message("Skipped line: " + line);
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
                        double value = std::stod(value_str);
                        if (key == "dt") dt = value;
                        params[key] = value;
                        log_message("Parsed parameter: " + key + "=" + std::to_string(value));
                    } catch (...) {
                        log_message("ERROR: Invalid parameter value: " + token);
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
                eq = new CustomEquation(name, equation, params, false);
            } else if (type == "particle") {
                eq = new CustomEquation(name, equation, params, true);
            } else if (type == "original") {
                eq = new OriginalEquation(name, params);
            } else if (type == "refined") {
                eq = new RefinedEquation(name, params);
            } else if (type == "tensor") {
                eq = new TensorEquation(name, params);
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