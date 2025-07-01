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

// Thread-safe logger
class Logger {
private:
    std::mutex mutex;
    std::vector<std::string> buffer;
    const size_t max_buffer_size = 100;

public:
    void log(const std::string& message);
    void flush();
};

// Tensor implementation
class Tensor {
public:
    Tensor();
    void set_diagonal(long double value);
    long double data[3][3];
};

// Particle implementation
class Particle {
public:
    Particle();
    long double x, y, z;
    long double vx, vy, vz;
    long double ax, ay, az;
};

class Simulation;

class Equation {
public:
    virtual ~Equation() = default;
    virtual void compute(Simulation& sim) = 0;
    virtual std::string output(const Simulation& sim) const = 0;
    virtual std::string name() const = 0;
    virtual bool is_valid() const = 0;
    virtual void setVisScale(long double s) = 0;
    virtual long double getVisScale() const = 0;
    virtual void setVisColorIntensity(long double c) = 0;
    virtual long double getVisColorIntensity() const = 0;
};

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
    void run(int steps);
    void setVisScale(const std::string& eq_name, long double scale);
    void setVisColorIntensity(const std::string& eq_name, long double intensity);
    void reset();

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
};

class ExpressionParser {
public:
    ExpressionParser();
    long double evaluate(const std::string& expr, const std::map<std::string, long double>& vars);
    bool validate(const std::string& expr) const;

private:
    std::map<std::string, std::function<long double(long double)>> functions;
    std::map<std::string, std::function<long double(long double, long double)>> binary_functions;
    std::map<std::string, long double> vars;
    std::string expression;
    size_t pos;

    std::string preprocess(const std::string& expr);
    void skipWhitespace();
    long double parseExpression();
    long double parseTerm();
    long double parseFactor();
    long double parseNumber();
    std::string parseToken();
};

class CustomEquation : public Equation {
public:
    CustomEquation(const std::string& name, const std::string& eq, const std::map<std::string, long double>& params,
                   const std::map<std::string, std::string>& derived_params, bool is_particle_based);
    bool is_valid() const override;
    void compute(Simulation& sim) override;
    std::string output(const Simulation& sim) const override;
    std::string name() const override;
    void setVisScale(long double s) override;
    long double getVisScale() const override;
    void setVisColorIntensity(long double c) override;
    long double getVisColorIntensity() const override;

private:
    std::string eq_name;
    std::string equation;
    ExpressionParser parser;
    std::map<std::string, long double> parameters;
    std::map<std::string, std::string> derived_parameters;
    bool is_particle_based;
    bool valid;
    long double scale, mass, c, G, mass1, mass2, r;
    long double wave_speed, freq, h, planck_freq, a, T, M, radius;
    long double l, m, n, l_quant, m_quant, gw_amp, gw_freq;
    long double E, B, kg_mass, spin, viscosity, temp, dim, alpha;
    long double H0, H1D, D_t, eps_5D, I0, D0, FiveD0, k, w, a0;
    long double vis_scale, vis_color_intensity;
    const long double max_ld;
    const long double min_ld;
    const long double duration;
    std::chrono::steady_clock::time_point start_time;
};

class TensorEquation : public Equation {
public:
    TensorEquation(const std::string& name, const std::string& eq, const std::map<std::string, long double>& params,
                   const std::map<std::string, std::string>& derived_params);
    bool is_valid() const override;
    void compute(Simulation& sim) override;
    std::string output(const Simulation& sim) const override;
    std::string name() const override;
    void setVisScale(long double s) override;
    long double getVisScale() const override;
    void setVisColorIntensity(long double c) override;
    long double getVisColorIntensity() const override;

private:
    std::string eq_name;
    std::string equation;
    std::map<std::string, long double> parameters;
    std::map<std::string, std::string> derived_parameters;
    long double H1D, D_t_coeff, eps_5D, alpha, I0, D0, FiveD0;
    long double vis_scale, vis_color_intensity;
    const long double max_ld;
    const long double min_ld;
    const long double duration;
    std::chrono::steady_clock::time_point start_time;
};

class QuantumEquation : public Equation {
public:
    QuantumEquation(const std::string& name, const std::map<std::string, long double>& params);
    bool is_valid() const override;
    void compute(Simulation& sim) override;
    std::string output(const Simulation& sim) const override;
    std::string name() const override;
    void setVisScale(long double s) override;
    long double getVisScale() const override;
    void setVisColorIntensity(long double c) override;
    long double getVisColorIntensity() const override;

private:
    std::string eq_name;
    long double hbar, m, a0;
    long double vis_scale, vis_color_intensity;
};

class OriginalEquation : public Equation {
public:
    OriginalEquation(const std::string& name, const std::map<std::string, long double>& params);
    bool is_valid() const override;
    void compute(Simulation& sim) override;
    std::string output(const Simulation& sim) const override;
    std::string name() const override;
    void setVisScale(long double s) override;
    long double getVisScale() const override;
    void setVisColorIntensity(long double c) override;
    long double getVisColorIntensity() const override;

private:
    std::string eq_name;
    long double H1D, D_t_coeff, eps_5D;
    long double vis_scale, vis_color_intensity;
};

class RefinedEquation : public Equation {
public:
    RefinedEquation(const std::string& name, const std::map<std::string, long double>& params);
    bool is_valid() const override;
    void compute(Simulation& sim) override;
    std::string output(const Simulation& sim) const override;
    std::string name() const override;
    void setVisScale(long double s) override;
    long double getVisScale() const override;
    void setVisColorIntensity(long double c) override;
    long double getVisColorIntensity() const override;

private:
    std::string eq_name;
    long double phi_base, G, k4, k5;
    long double vis_scale, vis_color_intensity;
    const long double max_ld;
    const long double min_ld;
    const long double duration;
    std::chrono::steady_clock::time_point start_time;
};

#endif