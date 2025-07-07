#include "equations.h"

// Global logger instance
Logger logger;

bool CustomEquation::is_valid() const {
    return valid;
}

std::string CustomEquation::output(const Simulation& sim) const {
    std::stringstream ss;
    ss << "# CustomEquation: " << eq_name << "\n";
    ss << "# t = " << sim.t << "\n";
    if (is_particle_based) {
        for (size_t i = 0; i < sim.particles.size(); ++i) {
            const auto& p = sim.particles[i];
            ss << "Particle " << i << ": x=" << p.x << ", y=" << p.y << ", z=" << p.z
               << ", vx=" << p.vx << ", vy=" << p.vy << ", vz=" << p.vz << "\n";
        }
    } else {
        for (int i = 0; i < sim.size; ++i) {
            for (int j = 0; j < sim.size; ++j) {
                for (int k = 0; k < sim.size; ++k) {
                    ss << "Grid [" << i << "][" << j << "][" << k << "]: scalar=" << sim.computed_scalar[i][j][k];
                    if (compute_divergence) {
                        ss << ", divergence=" << sim.divergence[i][j][k].data[0][0];
                    }
                    ss << "\n";
                }
            }
        }
    }
    return ss.str();
}

std::string CustomEquation::name() const {
    return eq_name;
}

void CustomEquation::setVisScale(long double s) {
    vis_scale = s;
}

long double CustomEquation::getVisScale() const {
    return vis_scale;
}

void CustomEquation::setVisColorIntensity(long double c) {
    vis_color_intensity = c;
}

long double CustomEquation::getVisColorIntensity() const {
    return vis_color_intensity;
}

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
                                if (compute_divergence && i > 0 && i < sim.size - 1 && j > 0 && j < sim.size - 1 && k > 0 && k < sim.size - 1) {
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

bool QuantumEquation::is_valid() const {
    return hbar > 0.0L && m > 0.0L && a0 > 0.0L;
}

std::string QuantumEquation::output(const Simulation& sim) const {
    std::stringstream ss;
    ss << "# QuantumEquation: " << eq_name << "\n";
    ss << "# t = " << sim.t << "\n";
    for (int i = 0; i < sim.size; ++i) {
        for (int j = 0; j < sim.size; ++j) {
            for (int k = 0; k < sim.size; ++k) {
                ss << "Grid [" << i << "][" << j << "][" << k << "]: scalar=" << sim.computed_scalar[i][j][k] << "\n";
            }
        }
    }
    return ss.str();
}

std::string QuantumEquation::name() const {
    return eq_name;
}

void QuantumEquation::setVisScale(long double s) {
    vis_scale = s;
}

long double QuantumEquation::getVisScale() const {
    return vis_scale;
}

void QuantumEquation::setVisColorIntensity(long double c) {
    vis_color_intensity = c;
}

long double QuantumEquation::getVisColorIntensity() const {
    return vis_color_intensity;
}

void QuantumEquation::compute(Simulation& sim) {
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
                        sim.scalar_4d[i][j][k][w] = vis_value * (1.0L + 0.1L * std::sin(w_coord + sim.t));
                    }
                }
            }
        }
    }
    sim.t += sim.dt;
    logger.log("Computed quantum step at t=" + std::to_string(sim.t));
    logger.flush();
}

void Simulation::load_equations() {
    equations.clear();
    std::map<std::string, long double> params;
    std::map<std::string, std::string> derived_params;

    // Example CustomEquation for a gravitational field
    params["G"] = 6.67430e-11L;
    params["mass1"] = 1.989e30L;
    params["mass2"] = 5.972e24L;
    params["r"] = 1.496e11L;
    equations.emplace_back(std::make_unique<CustomEquation>(
        "Gravitational", "-G*mass1*mass2/(r*r)", params, derived_params, true));

    // Example CustomEquation for a tensor field with divergence
    params.clear();
    params["H"] = 1.0e-18L;
    params["eps"] = 0.1L;
    params["L0"] = 1.0L;
    params["alpha"] = 0.01L;
    equations.emplace_back(std::make_unique<CustomEquation>(
        "TensorField", "sum(H*exp(-alpha*d*t)/(d*d*d*d))", params, derived_params, false, true));

    // Example CustomEquation for a refined field
    params.clear();
    params["phi"] = 1.0L;
    params["G"] = 1.0L;
    params["k4"] = 0.1L;
    params["k5"] = 0.01L;
    params["H"] = 1.0e-18L;
    params["eps"] = 0.1L;
    params["L0"] = 1.0L;
    equations.emplace_back(std::make_unique<CustomEquation>(
        "RefinedField", "phi*sum(H*exp(-k4*d*t)/(d*d*d*d))", params, derived_params, false));

    // QuantumEquation for hydrogen atom
    params.clear();
    params["hbar"] = 1.05457182e-34L;
    params["m"] = 9.1093837e-31L;
    params["a0"] = 5.29e-11L;
    equations.emplace_back(std::make_unique<QuantumEquation>("HydrogenAtom", params));
}

void Simulation::init(int grid_size, long double dx_val, long double dt_val, int num_particles) {
    size = grid_size;
    dx = dx_val;
    dt = dt_val;
    t = 0.0L;
    current_equation = 0;

    ricci.resize(size, std::vector<std::vector<Tensor>>(size, std::vector<Tensor>(size)));
    divergence.resize(size, std::vector<std::vector<Tensor>>(size, std::vector<Tensor>(size)));
    computed_scalar.resize(size, std::vector<std::vector<long double>>(size, std::vector<long double>(size, 0.0L)));
    scalar_4d.resize(size, std::vector<std::vector<std::vector<long double>>>(
        size, std::vector<std::vector<long double>>(size, std::vector<long double>(size, 0.0L))));

    particles.clear();
    particles.resize(num_particles);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<long double> dist(-size * dx / 2, size * dx / 2);
    for (auto& p : particles) {
        p.x = dist(gen);
        p.y = dist(gen);
        p.z = dist(gen);
        p.vx = dist(gen) * 0.1L;
        p.vy = dist(gen) * 0.1L;
        p.vz = dist(gen) * 0.1L;
    }
}

void Simulation::run(int steps, const std::string& eq_name) {
    for (size_t i = 0; i < equations.size(); ++i) {
        if (equations[i]->name() == eq_name) {
            current_equation = i;
            break;
        }
    }
    for (int i = 0; i < steps; ++i) {
        step();
    }
}

std::string Simulation::output(const std::string& eq_name) const {
    for (const auto& eq : equations) {
        if (eq->name() == eq_name) {
            return eq->output(*this);
        }
    }
    return "Equation not found: " + eq_name;
}

void Simulation::setVisParameters(const std::string& eq_name, long double scale, long double color_intensity) {
    for (const auto& eq : equations) {
        if (eq->name() == eq_name) {
            eq->setVisScale(scale);
            eq->setVisColorIntensity(color_intensity);
        }
    }
}

void Simulation::setVisScale(const std::string& eq_name, long double scale) {
    for (const auto& eq : equations) {
        if (eq->name() == eq_name) {
            eq->setVisScale(scale);
        }
    }
}

void Simulation::setVisColorIntensity(const std::string& eq_name, long double intensity) {
    for (const auto& eq : equations) {
        if (eq->name() == eq_name) {
            eq->setVisColorIntensity(intensity);
        }
    }
}

void Simulation::reset() {
    t = 0.0L;
    for (auto& vec : computed_scalar) {
        for (auto& row : vec) {
            std::fill(row.begin(), row.end(), 0.0L);
        }
    }
    for (auto& vec : ricci) {
        for (auto& row : vec) {
            for (auto& tensor : row) {
                std::fill(&tensor.data[0][0], &tensor.data[0][0] + 9, 0.0L);
            }
        }
    }
    for (auto& vec : divergence) {
        for (auto& row : vec) {
            for (auto& tensor : row) {
                std::fill(&tensor.data[0][0], &tensor.data[0][0] + 9, 0.0L);
            }
        }
    }
    for (auto& vec : scalar_4d) {
        for (auto& row : vec) {
            for (auto& col : row) {
                std::fill(col.begin(), col.end(), 0.0L);
            }
        }
    }
    for (auto& p : particles) {
        p.x = p.y = p.z = p.vx = p.vy = p.vz = p.ax = p.ay = p.az = 0.0L;
    }
}