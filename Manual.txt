Overview
This physics visualization tool allows physicists and mathematicians to simulate and visualize a variety of physical systems, ranging from classical mechanics (e.g., projectile motion, harmonic oscillators) to general relativity (e.g., Ricci tensors, Schwarzschild metric) and quantum mechanics (e.g., hydrogen wavefunction). The tool supports both particle-based simulations and grid-based scalar field visualizations, controlled through a configuration file (equations.txt) and an interactive 3D OpenGL interface.
Straightforward Usage Instructions
1. Installation and Setup
Prerequisites:
# Dependencies to install manually (run once):
sudo apt update
sudo apt install -y make cmake build-essential libsdl2-dev libglu1-mesa-dev mingw-w64 wget unzip libsdl2-ttf-dev

The tool is compatible with Linux, Android, and Windows (with appropriate SDL2 setup).

Running the Tool:
Place your equations.txt file in the same directory as the executable.
Ensure you have an output folder for the calculations
Run the program: ./simulation or simulation.exe or simulation on Android.

Physics Visualization Controls:

- Esc: Quit the program
- 1: Switch to POINTS display mode (shows scalar field as colored points)
- 2: Switch to ISOSURFACE display mode (shows scalar field as 3D surfaces)
- 3: Switch to WIREFRAME display mode (shows scalar field as a wireframe grid)
- 4: Switch to PARTICLES display mode (shows particle positions and motion)
- 5: Switch to HYBRID display mode (shows both scalar field and particles)
- Left Arrow: Switch to the previous equation
- Right Arrow: Switch to the next equation
- R: Reset simulation (time to 0, camera to default position and zoom)
- P: Clear particles and reset simulation
- Space: Toggle camera rotation (pause/resume automatic rotation)
- Up Arrow: Zoom in (move camera closer to the simulation)
- Down Arrow: Zoom out (move camera farther from the simulation)
- Z: Increase grid size (higher resolution, max 50)
- X: Decrease grid size (lower resolution, min 5)
- C: Center camera (reset angle and zoom to default)

2. Configuring Equations (equations.txt)
File Structure:
The tool reads equations and parameters from equations.txt.

Lines starting with # are comments and ignored.

The file supports three types of entries:
Simulation Parameters:
Format: param key=value key=value ...

Example: param dt=0.016 g=9.8

Supported parameters:
dt: Time step for the simulation (default: 0.016 seconds).

Other parameters (e.g., g) can be used in equations.

Particle-Based Equations:
Format:
equation <name>
<variable> = <expression>
<variable> = <expression>
...

Supported variables: x, y, z (position), vx, vy, vz (velocity), ax, ay, az (acceleration), t (time), dt.

Supported operations: +, -, *, /, sin, cos, sqrt, exp.

Example:
equation Projectile
x = x + vx * dt
y = y + vy * dt
vx = vx
vy = vy - g * dt

Grid-Based Equations:
Format: <equation> <type> [parameters]

Types:
Original Equations: For tensor-based computations (e.g., Ricci tensors).
Identified by the presence of ∇ or \nabla.

Parameters: H1D, D_t, eps_5D.

Example: ∇μR3Dμν(t)=H1D(α)gμν+I3Dμν−4D(t)+ϵ5Dμν(t) H1D=1.0 D_t=0.1 eps_5D=0.01

Refined Equations: For curvature-like computations.
Format: <name> refined phi=<value> G=<value> k4=<value> k5=<value>

Example: Refined_Eq1 refined phi=1.0 G=1.0 k4=0.1 k5=0.01

Custom Equations: Predefined physics equations.
Format: <name> custom [parameters]

Supported equations and parameters:
E_equals_mc2: mass, c

Newton_Gravity: G, mass1, mass2, r

Wave_Equation: wave_speed, freq

Planck_Energy: h, freq

Kepler_Third: a, T

Schwarzschild: M, r

Spherical_Harmonic: l, m

Hydrogen_Wavefunction: n, l, m

Gravitational_Wave: amp, freq

Electromagnetic_Field: E, B

Klein_Gordon: mass

Dirac_Field: spin

Navier_Stokes: viscosity

Blackbody_Radiation: T

String_Theory: dim, alpha

Cosmological_Expansion: H0

Example: E_equals_mc2 custom mass=1.0 c=299792458

3. Running the Visualization
Window:
A window titled "Physics Visualization" (800x600) opens, displaying the simulation.

The current equation name is shown in the top-left corner.

4. Output
Simulation Data:
Output files are generated in the output/ directory as output/equation_<name>.txt.

For particle-based equations: Lists particle positions and velocities.

For grid-based equations: Lists scalar field values at each grid point.

Logs:
Program logs are written to log.txt, including initialization, equation switches, and errors.

Deep-Dive for Advanced Users
Simulation Mechanics
Hybrid Simulation:
The tool supports two simulation types:
Particle-Based:
Used for equations defined with equation blocks (e.g., Projectile).

Simulates 10 particles by default, each with position (x, y, z), velocity (vx, vy, vz), and acceleration (ax, ay, az).

Equations are evaluated using a simple expression parser supporting basic math operations and functions (sin, cos, sqrt, exp).

Particles are updated using Euler integration: x += vx * dt + 0.5 * ax * dt^2, vx += ax * dt, etc.

Grid-Based:
Used for Original, Refined, and Custom equations (e.g., E_equals_mc2).

Simulates a 10x10x10 grid with spacing dx=1.0 (configurable during initialization).

Computes a 3D scalar field (scalar) and a 4D scalar field (scalar_4d) with a w dimension.

OriginalEquation: Computes tensor divergences with time-dependent terms.

RefinedEquation: Computes Ricci tensors with gravitational and perturbation terms.

CustomEquation: Computes predefined physics equations (e.g., E=mc^2, Newton’s gravity) as scalar values across the grid.

Time Evolution:
The simulation updates every 16 ms (60 FPS).

Time step dt is set via equations.txt (default: 0.016 s).

Time t increments by dt each frame and is reset to 0 when switching equations.

Visualization Details
Grid-Based Visualization:
Points Mode:
Each grid point is a dot, colored by:
Red: Normalized scalar value.

Green: Normalized scalar_4d value at the current w slice.

Blue: 1 - red.

The w slice animates over time: w_slice = (t * 2) % grid_size.

Isosurface Mode:
Uses Marching Cubes to render an isosurface at the midpoint of the scalar range.

Colors vary based on scalar_4d values.

Includes lighting for depth.

Wireframe Mode:
Draws lines between grid points, offset by scalar_4d values for a dynamic effect.

Particle-Based Visualization:
Particles are rendered as points, colored by velocity magnitude (red increases with speed, blue decreases).

Hybrid Mode:
Combines particle rendering with grid points.

Camera:
Orbits at a radius of 20 units, with camera_angle incrementing by 0.01 radians per frame.

Axes:
X (red), Y (green), Z (blue) axes extend from (0,0,0) to (10,0,0), (0,10,0), and (0,0,10), respectively.

Limitations and Considerations
Grid Size:
Fixed at 10x10x10 by default. Increasing this requires modifying main.cpp (sim.initialize).

Expression Parser:
The parser for particle-based equations is basic and may fail on complex expressions. Stick to simple operations.

Performance:
Grid-based equations can be computationally intensive for large grids.

The 60 FPS update rate may slow down on lower-end hardware.

Units:
Assumes SI units (meters, seconds) for particle-based equations. Grid-based equations are abstract and may not correspond to physical units.

Extending the Tool
Adding New Custom Equations:
Modify CustomEquation::compute in equations.cpp to add new predefined equations.

Customizing Visualization:
Adjust renderer.cpp to change colors, camera behavior, or add new display modes.

Improving Performance:
Optimize the Marching Cubes algorithm in ISOSURFACE mode or reduce the grid size for faster rendering.
