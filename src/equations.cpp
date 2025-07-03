#include "equations.h"

// Global logger instance
Logger logger;

void CustomEquation::compute(Simulation& sim) {
    if (!valid) {
        logger.log("ERROR: Cannot compute invalid equation: " + eq_name);
        sim.t += sim.dt;
        return;
    }
    auto current_time = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::duration<long double>>(current_time - start_time).count();

    if (is_particle_based) {
        #pragma omp parallel for
        for (size_t i = 0; i < sim.particles.size(); ++i) {
            std::map<std::string, long double> vars = parameters;
            vars["x"] = sim.particles[i].x;
            vars["y"] = sim.particles[i].y;
            vars["z"] = sim.particles[i].z;
            vars["vx"] = sim.particles[i].vx;
            vars["vy"] = sim.particles[i].vy;
            vars["vz"] = sim.particles[i].vz;
            vars["ax"] = sim.particles[i].ax;
            vars["ay"] = sim.particles[i].ay;
            vars["az"] = sim.particles[i].az;
            vars["t"] = sim.t;
            vars["dt"] = sim.dt;
            vars["r"] = std::sqrt(sim.particles[i].x * sim.particles[i].x + sim.particles[i].y * sim.particles[i].y + sim.particles[i].z * sim.particles[i].z);
            vars["theta"] = vars["r"] > 1e-6L ? std::acos(sim.particles[i].z / vars["r"]) : 0.0L;
            vars["phi"] = std::atan2(sim.particles[i].y, sim.particles[i].x);
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

            ExpressionParser local_parser;
            for (const auto& dp : derived_parameters) {
                try {
                    vars[dp.first] = local_parser.evaluate(dp.second, vars);
                } catch (const std::exception& e) {
                    logger.log("ERROR: Failed to evaluate derived parameter '" + dp.first + "': " + e.what());
                    vars[dp.first] = 0.0L;
                }
            }

            try {
                long double value = equation.find("sum(") == 0 ? compute_sum_field(elapsed, vars, equation) : local_parser.evaluate(equation, vars);
                long double vis_value = value * scale * vis_scale;
                vis_value = std::tanh(vis_value / vis_color_intensity);
                long double ax_new = vis_value;
                long double ay_new = vis_value;
                long double az_new = vis_value;

                #pragma omp critical
                {
                    sim.particles[i].x += sim.particles[i].vx * sim.dt + 0.5L * ax_new * sim.dt * sim.dt;
                    sim.particles[i].y += sim.particles[i].vy * sim.dt + 0.5L * ay_new * sim.dt * sim.dt;
                    sim.particles[i].z += sim.particles[i].vz * sim.dt + 0.5L * az_new * sim.dt * sim.dt;
                    sim.particles[i].vx += ax_new * sim.dt;
                    sim.particles[i].vy += ay_new * sim.dt;
                    sim.particles[i].vz += az_new * sim.dt;
                    sim.particles[i].ax = ax_new;
                    sim.particles[i].ay = ay_new;
                    sim.particles[i].az = az_new;
                }
            } catch (const std::exception& e) {
                logger.log("ERROR: Failed to evaluate equation '" + equation + "': " + e.what());
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
                        long double x = static_cast<long double>(i - sim.size / 2) * sim.dx;
                        long double y = static_cast<long double>(j - sim.size / 2) * sim.dx;
                        long double z = static_cast<long double>(k - sim.size / 2) * sim.dx;
                        long double r = std::sqrt(x * x + y * y + z * z);

                        std::map<std::string, long double> vars = parameters;
                        vars["x"] = x;
                        vars["y"] = y;
                        vars["z"] = z;
                        vars["r"] = r;
                        vars["t"] = sim.t;
                        vars["theta"] = r > 1e-6L ? std::acos(z / r) : 0.0L;
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
                                logger.log("ERROR: Failed to evaluate derived parameter '" + dp.first + "': " + e.what());
                                vars[dp.first] = 0.0L;
                            }
                        }

                        try {
                            long double value = equation.find("sum(") == 0 ? compute_sum_field(elapsed, vars, equation) : local_parser.evaluate(equation, vars);
                            long double vis_value = value * scale * vis_scale;
                            vis_value = std::tanh(vis_value / vis_color_intensity);
                            #pragma omp critical
                            {
                                sim.ricci[i][j][k].set_diagonal(vis_value);
                                sim.computed_scalar[i][j][k] = vis_value;
                                for (int w = 0; w < sim.size; ++w) {
                                    long double w_coord = static_cast<long double>(w - sim.size / 2) * sim.dx;
                                    sim.scalar_4d[i][j][k][w] = vis_value * (1.0L + 0.3L * std::sin(w_coord + sim.t));
                                }
                            }
                        } catch (const std::exception& e) {
                            logger.log("ERROR: Failed to evaluate equation '" + equation + "': " + e.what());
                        }
                    }
                }
            }
        }
    }
    sim.t += sim.dt;
    logger.log("Computed step at t=" + std::to_string(sim.t) + ", progress=" + std::to_string(std::min(elapsed / duration, 1.0L)));
    logger.flush();
}

void TensorEquation::compute(Simulation& sim) {
    auto current_time = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::duration<long double>>(current_time - start_time).count();

    #pragma omp parallel
    {
        ExpressionParser local_parser;
        #pragma omp for collapse(3)
        for (int i = 0; i < sim.size; ++i) {
            for (int j = 0; j < sim.size; ++j) {
                for (int k = 0; k < sim.size; ++k) {
                    long double x = static_cast<long double>(i - sim.size / 2) * sim.dx;
                    long double y = static_cast<long double>(j - sim.size / 2) * sim.dx;
                    long double z = static_cast<long double>(k - sim.size / 2) * sim.dx;
                    long double r = std::sqrt(x * x + y * y + z * z);

                    std::map<std::string, long double> vars = parameters;
                    vars["x"] = x;
                    vars["y"] = y;
                    vars["z"] = z;
                    vars["r"] = r;
                    vars["t"] = sim.t;
                    vars["theta"] = r > 1e-6L ? std::acos(z / r) : 0.0L;
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
                            logger.log("ERROR: Failed to evaluate derived parameter '" + dp.first + "': " + e.what());
                            vars[dp.first] = 0.0L;
                        }
                    }

                    try {
                        long double value = equation.find("sum(") == 0 ? compute_sum_field(elapsed, vars, equation) : local_parser.evaluate(equation, vars);
                        long double vis_value = value * vis_scale;
                        vis_value = std::tanh(vis_value / vis_color_intensity);
                        #pragma omp critical
                        {
                            sim.ricci[i][j][k].set_diagonal(vis_value);
                            sim.computed_scalar[i][j][k] = vis_value;

                            if (i > 0 && i < sim.size - 1 && j > 0 && j < sim.size - 1 && k > 0 && k < sim.size - 1) {
                                Tensor& div = sim.divergence[i][j][k];
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
                            }

                            for (int w = 0; w < sim.size; ++w) {
                                long double w_coord = static_cast<long double>(w - sim.size / 2) * sim.dx;
                                sim.scalar_4d[i][j][k][w] = vis_value * (1.0L + 0.1L * std::sin(w_coord + sim.t));
                            }
                        }
                    } catch (const std::exception& e) {
                        logger.log("ERROR: Failed to evaluate equation '" + equation + "': " + e.what());
                    }
                }
            }
        }
    }
    sim.t += sim.dt;
    logger.log("Computed step at t=" + std::to_string(sim.t) + ", progress=" + std::to_string(std::min(elapsed / duration, 1.0L)));
    logger.flush();
}

void RefinedEquation::compute(Simulation& sim) {
    auto current_time = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::duration<long double>>(current_time - start_time).count();

    #pragma omp parallel for collapse(3)
    for (int i = 0; i < sim.size; ++i) {
        for (int j = 0; j < sim.size; ++j) {
            for (int k = 0; k < sim.size; ++k) {
                long double x = static_cast<long double>(i - sim.size / 2) * sim.dx;
                long double y = static_cast<long double>(j - sim.size / 2) * sim.dx;
                long double z = static_cast<long double>(k - sim.size / 2) * sim.dx;
                long double r = std::sqrt(x * x + y * y + z * z) + 1e-6L;

                std::map<std::string, long double> vars;
                vars["x"] = x;
                vars["y"] = y;
                vars["z"] = z;
                vars["r"] = r;
                vars["t"] = sim.t;
                vars["phi"] = phi_base;
                vars["G"] = G;
                vars["k4"] = k4;
                vars["k5"] = k5;
                vars["H"] = 1.0e-18L;
                vars["eps"] = 0.1L;
                vars["L0"] = 1.0L;
                vars["alpha"] = 0.01L;

                long double field_value = compute_sum_field(elapsed, vars, "");
                long double phi = phi_base * (1.0L + k4 * field_value);
                long double vis_value = phi * vis_scale;
                vis_value = std::tanh(vis_value / vis_color_intensity);

                #pragma omp critical
                {
                    sim.computed_scalar[i][j][k] = vis_value;
                    sim.ricci[i][j][k].set_diagonal(vis_value);
                    for (int w = 0; w < sim.size; ++w) {
                        long double w_coord = static_cast<long double>(w - sim.size / 2) * sim.dx;
                        sim.scalar_4d[i][j][k][w] = vis_value * (1.0L + k5 * std::sin(w_coord + sim.t));
                    }
                }
            }
        }
    }
    sim.t += sim.dt;
    logger.log("Computed step at t=" + std::to_string(sim.t) + ", progress=" + std::to_string(std::min(elapsed / duration, 1.0L)));
    logger.flush();
}

void Simulation::init(int grid_size, long double dx_val, long double dt_val, int num_particles) {
    size = grid_size;
    dx = dx_val;
    dt = dt_val;
    t = 0.0L;
    current_equation = 0;

    particles.clear();
    particles.resize(num_particles);
    std::mt19937 gen(static_cast<unsigned>(std::time(nullptr)));
    std::uniform_real_distribution<long double> dist(-1.0L, 1.0L);
    for (auto& p : particles) {
        p.x = dist(gen) * size * dx / 2;
        p.y = dist(gen) * size * dx / 2;
        p.z = dist(gen) * size * dx / 2;
        p.vx = dist(gen) * 0.1L;
        p.vy = dist(gen) * 0.1L;
        p.vz = dist(gen) * 0.1L;
    }

    ricci.assign(size, std::vector<std::vector<Tensor>>(size, std::vector<Tensor>(size)));
    divergence.assign(size, std::vector<std::vector<Tensor>>(size, std::vector<Tensor>(size)));
    computed_scalar.assign(size, std::vector<std::vector<long double>>(size, std::vector<long double>(size, 0.0L)));
    scalar_4d.assign(size, std::vector<std::vector<std::vector<long double>>>(
        size, std::vector<std::vector<long double>>(size, std::vector<long double>(size, 0.0L))));
}

void Simulation::load_equations() {
    equations.clear();
    std::map<std::string, long double> params;
    std::map<std::string, std::string> derived_params;

    params["H1D"] = 1.0L;
    params["D_t"] = 0.1L;
    params["eps_5D"] = 0.1L;
    equations.push_back(std::make_unique<OriginalEquation>("Original", params));

    params.clear();
    params["hbar"] = 1.05457182e-34L;
    params["m"] = 9.1093837e-31L;
    params["a0"] = 5.29e-11L;
    equations.push_back(std::make_unique<QuantumEquation>("Quantum", params));

    params.clear();
    params["scale"] = 1.0L;
    params["c"] = 299792458.0L;
    equations.push_back(std::make_unique<CustomEquation>("Wave", "sin(2*pi*freq*t)", params, derived_params, false));

    params.clear();
    params["G"] = 6.67430e-11L;
    params["mass1"] = 1.989e30L;
    params["mass2"] = 1.989e30L;
    params["r"] = 1.0e6L;
    derived_params["F"] = "G*mass1*mass2/(r*r)";
    equations.push_back(std::make_unique<CustomEquation>("Gravitational Force", "F/r", params, derived_params, true));

    params.clear();
    params["H"] = 1.0e-18L;
    params["eps"] = 0.1L;
    params["L0"] = 1.0L;
    params["alpha"] = 0.01L;
    equations.push_back(std::make_unique<TensorEquation>("SumField", "sum(H,eps,L0)", params, derived_params));

    params.clear();
    params["phi"] = 1.0L;
    params["G"] = 1.0L;
    params["k4"] = 0.1L;
    params["k5"] = 0.01L;
    equations.push_back(std::make_unique<RefinedEquation>("Refined", params));

    logger.log("Loaded " + std::to_string(equations.size()) + " equations");
}

void Simulation::run(int steps, const std::string& eq_name) {
    size_t target_eq = current_equation;
    if (!eq_name.empty()) {
        for (size_t i = 0; i < equations.size(); ++i) {
            if (equations[i]->name() == eq_name) {
                target_eq = i;
                break;
            }
        }
    }
    if (target_eq < equations.size()) {
        current_equation = target_eq;
        for (int i = 0; i < steps; ++i) {
            equations[current_equation]->compute(*this);
            logger.log("Step " + std::to_string(i + 1) + "/" + std::to_string(steps) + " at t=" + std::to_string(t));
        }
    } else {
        logger.log("ERROR: Equation '" + eq_name + "' not found");
    }
    logger.flush();
}

std::string Simulation::output(const std::string& eq_name) const {
    if (equations.empty()) return "# No equations loaded\n";
    size_t target_eq = current_equation;
    if (!eq_name.empty()) {
        for (size_t i = 0; i < equations.size(); ++i) {
            if (equations[i]->name() == eq_name) {
                target_eq = i;
                break;
            }
        }
    }
    if (target_eq < equations.size()) {
        return equations[target_eq]->output(*this);
    }
    return "# Equation '" + eq_name + "' not found\n";
}

void Simulation::setVisParameters(const std::string& eq_name, long double scale, long double color_intensity) {
    for (auto& eq : equations) {
        if (eq_name.empty() || eq->name() == eq_name) {
            eq->setVisScale(scale);
            eq->setVisColorIntensity(color_intensity);
            logger.log("Set visualization parameters for " + eq->name() + ": scale=" + std::to_string(scale) +
                       ", color_intensity=" + std::to_string(color_intensity));
        }
    }
    logger.flush();
}

void Simulation::setVisScale(const std::string& eq_name, long double scale) {
    for (auto& eq : equations) {
        if (eq_name.empty() || eq->name() == eq_name) {
            eq->setVisScale(scale);
            logger.log("Set visualization scale for " + eq->name() + ": scale=" + std::to_string(scale));
        }
    }
    logger.flush();
}

void Simulation::setVisColorIntensity(const std::string& eq_name, long double intensity) {
    for (auto& eq : equations) {
        if (eq_name.empty() || eq->name() == eq_name) {
            eq->setVisColorIntensity(intensity);
            logger.log("Set visualization color intensity for " + eq->name() + ": intensity=" + std::to_string(intensity));
        }
    }
    logger.flush();
}

void Simulation::reset() {
    t = 0.0L;
    particles.clear();
    ricci.assign(size, std::vector<std::vector<Tensor>>(size, std::vector<Tensor>(size)));
    divergence.assign(size, std::vector<std::vector<Tensor>>(size, std::vector<Tensor>(size)));
    computed_scalar.assign(size, std::vector<std::vector<long double>>(size, std::vector<long double>(size, 0.0L)));
    scalar_4d.assign(size, std::vector<std::vector<std::vector<long double>>>(
        size, std::vector<std::vector<long double>>(size, std::vector<long double>(size, 0.0L))));
    logger.log("Simulation reset");
    logger.flush();
}