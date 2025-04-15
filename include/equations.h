#ifndef EQUATIONS_H
#define EQUATIONS_H

#include <string>
#include <vector>
#include <set>
#include <map>

// Declare log_message function
void log_message(const std::string& message);

// Tensor class for grid-based computations
class Tensor {
public:
    double data[3][3];
    Tensor();
    void set_diagonal(double value);
};

// Particle structure for particle-based simulations
struct Particle {
    double x, y, z;    // Position
    double vx, vy, vz; // Velocity
    double ax, ay, az; // Acceleration
    Particle();
};

// Abstract base class for equations
class Equation {
public:
    virtual ~Equation() {}
    virtual void compute(class Simulation& sim) = 0;
    virtual std::string output(const Simulation& sim) const = 0;
    virtual std::string name() const = 0;
};

// Forward declarations
class OriginalEquation;
class RefinedEquation;
class CustomEquation;

// Simulation class to manage both particles and grid
class Simulation {
public:
    // Particle-based data
    std::vector<Particle> particles;
    // Grid-based data
    int size;
    double dx;
    std::vector<std::vector<std::vector<Tensor>>> ricci;
    std::vector<std::vector<std::vector<Tensor>>> divergence;
    std::vector<std::vector<std::vector<double>>> scalar;
    std::vector<std::vector<std::vector<double>>> computed_scalar;
    std::vector<std::vector<std::vector<std::vector<double>>>> scalar_4d;

    std::vector<Equation*> equations;
    double t;
    double dt;
    int current_equation;
    std::set<std::string> output_written;

    Simulation();
    void initialize(int num_particles, int grid_size, double grid_dx, double dt_val);
    void load_equations(const std::string& filename);
    void compute();
    void save_output();
    void switch_equation(int index);

    // Added methods for main.cpp compatibility
    void setup();
    void step();
    void prev_equation();
    void next_equation();
};

#endif