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
    ExpressionParser();
    double evaluate(const std::string& expr, const std::map<std::string, double>& vars);
    bool validate(const std::string& expr) const;

private:
    std::map<std::string, std::function<double(double)>> functions;
    std::map<std::string, std::function<double(double, double)>> binary_functions;
    std::map<std::string, double> vars;
    std::string expression;
    size_t pos;

    std::string preprocess(const std::string& expr);
    void skipWhitespace();
    double parseExpression();
    double parseTerm();
    double parseFactor();
    double parseNumber();
    std::string parseToken();
};

class CustomEquation : public Equation {
public:
    CustomEquation(const std::string& name, const std::string& eq, const std::map<std::string, double>& params,
                   const std::map<std::string, std::string>& derived_params, bool is_particle_based);
    bool is_valid() const override;
    void compute(Simulation& sim) override;
    std::string output(const Simulation& sim) const override;
    std::string name() const override;
    void setVisScale(double s);
    double getVisScale() const;
    void setVisColorIntensity(double c);
    double getVisColorIntensity() const;

private:
    std::string eq_name;
    std::string equation;
    ExpressionParser parser;
    std::map<std::string, double> parameters;
    std::map<std::string, std::string> derived_parameters;
    bool is_particle_based;
    bool valid;
    double scale, mass, c, G, mass1, mass2, r, wave_speed, freq, h, planck_freq, a, T, M, radius;
    double l, m, n, l_quant, m_quant, gw_amp, gw_freq, E, B, kg_mass, spin, viscosity, temp, dim, alpha, H0;
    double H1D, D_t, eps_5D, I0, D0, FiveD0, k, w, a0; // Added for tensor and other equations
    double vis_scale, vis_color_intensity;
};

class TensorEquation : public Equation {
public:
    TensorEquation(const std::string& name, const std::string& eq, const std::map<std::string, double>& params,
                   const std::map<std::string, std::string>& derived_params);
    void compute(Simulation& sim) override;
    std::string output(const Simulation& sim) const override;
    std::string name() const override;

private:
    std::string eq_name;
    std::string equation;
    ExpressionParser parser;
    std::map<std::string, double> parameters;
    std::map<std::string, std::string> derived_parameters;
    double H1D, D_t_coeff, eps_5D, alpha, I0, D0, FiveD0;
};

class QuantumEquation : public Equation {
public:
    QuantumEquation(const std::string& name, const std::map<std::string, double>& params);
    void compute(Simulation& sim) override;
    std::string output(const Simulation& sim) const override;
    std::string name() const override;

private:
    std::string eq_name;
    double hbar, m, a0;
};

class OriginalEquation : public Equation {
public:
    OriginalEquation(const std::string& name, const std::map<std::string, double>& params);
    void compute(Simulation& sim) override;
    std::string output(const Simulation& sim) const override;
    std::string name() const override;

private:
    std::string eq_name;
    double H1D, D_t_coeff, eps_5D;
};

class RefinedEquation : public Equation {
public:
    RefinedEquation(const std::string& name, const std::map<std::string, double>& params);
    void compute(Simulation& sim) override;
    std::string output(const Simulation& sim) const override;
    std::string name() const override;

private:
    std::string eq_name;
    double phi_base, G, k4, k5;
};

#endif