\# N-Link Pendulum Simulation



An interactive numerical simulation of a coupled \*\*N-link pendulum\*\*, developed to study nonlinear dynamics, coupled oscillations, energy transfer, and the numerical solution of the equations of motion.



The simulation derives and solves the coupled equations of motion for an arbitrary number of pendulum links and provides an interactive browser-based visualization of the resulting dynamics.



\---



\## Overview



A multi-link pendulum is a nonlinear dynamical system in which the motion of each link is coupled to every other link. Even though the governing equations are deterministic, the system can exhibit highly sensitive and complex motion as the number of links and initial energy increase.



This project implements an \*\*N-link pendulum solver\*\* using the Lagrangian formulation of mechanics and numerically integrates the resulting equations of motion using `SciPy`.



The simulation provides:



\* Real-time visualization of the pendulum

\* Adjustable number of links

\* Adjustable gravitational acceleration

\* Adjustable damping

\* Multiple initial configurations

\* Angle-versus-time plots

\* Kinetic, potential, and total energy plots

\* Animation controls and zooming

\* A built-in mathematical explanation of the model



\---



\## Physics



\### Generalized Coordinates



For an N-link pendulum, the angular displacement of each link is represented by



\[

\\theta\_1,\\theta\_2,\\ldots,\\theta\_N

]



where each angle is measured from the downward vertical.



The position of the end of link (i) is obtained from the cumulative contribution of all links above it:



\[

x\_i = \\sum\_{k=1}^{i} L\_k\\sin\\theta\_k

]



\[

y\_i = -\\sum\_{k=1}^{i} L\_k\\cos\\theta\_k

]



where:



\* (L\_k) = length of link (k)

\* (\\theta\_k) = angular displacement of link (k)



\---



\## Lagrangian Formulation



The equations of motion are obtained using the Lagrangian



\[

\\mathcal{L}=T-V

]



where:



\* (T) = total kinetic energy

\* (V) = total gravitational potential energy



For the system,



\[

T=\\frac{1}{2}\\sum\_{i=1}^{N}m\_i

\\left(\\dot{x}\_i^2+\\dot{y}\_i^2\\right)

]



and



\[

V=\\sum\_{i=1}^{N}m\_i g y\_i

]



The Euler-Lagrange equations are then used to obtain the coupled equations governing the angular coordinates.



\---



\## Matrix Form of the Equations of Motion



The resulting equations are expressed in the form



\[

M(\\theta)\\ddot{\\theta}=b(\\theta,\\dot{\\theta})

]



where:



\* (M(\\theta)) is the configuration-dependent mass matrix

\* (\\theta) is the vector of angular positions

\* (\\dot{\\theta}) is the vector of angular velocities

\* (b) contains the nonlinear, gravitational, and damping terms



For the implementation, the elements of the mass matrix are



\\left(

\\sum\_{i=\\max(j,k)}^{N}m\_i

\\right)

L\_jL\_k

\\cos(\\theta\_j-\\theta\_k)

]



The right-hand side contains the nonlinear coupling terms, gravitational terms, and damping:



\[

b\_j =

\-\\sum\_k

\\left(

\\sum\_{i=\\max(j,k)}^{N}m\_i

\\right)

L\_jL\_k

\\sin(\\theta\_j-\\theta\_k)\\dot{\\theta}\_k^2

]



\[

\-\\left(\\sum\_{i=j}^{N}m\_i\\right)gL\_j\\sin\\theta\_j

\-\\gamma\\dot{\\theta}\_j

]



The angular accelerations are obtained by solving



\[

\\ddot{\\theta}=M^{-1}b

]



rather than explicitly calculating the matrix inverse. The implementation uses a linear system solve:



```python

alpha = np.linalg.solve(M, b)

```



This gives the numerical angular accelerations required by the ODE solver.



\---



\## Numerical Method



The equations of motion form a system of coupled nonlinear ordinary differential equations.



The state vector is



\[

y =

\[

\\theta\_1,\\ldots,\\theta\_N,

\\dot{\\theta}\_1,\\ldots,\\dot{\\theta}\_N

]

]



The derivative of this state is



\[

\\dot{y} =

\[

\\dot{\\theta}\_1,\\ldots,\\dot{\\theta}\_N,

\\ddot{\\theta}\_1,\\ldots,\\ddot{\\theta}\_N

]

]



The system is numerically integrated using the `solve\_ivp` function from SciPy.



The current implementation uses the \*\*RK45\*\* method:



```python

solve\_ivp(

&#x20;   ...,

&#x20;   method="RK45",

&#x20;   rtol=1e-3,

&#x20;   atol=1e-3

)

```



RK45 is an adaptive Runge-Kutta method that adjusts its integration step according to the local numerical error.



\---



\## Damping



A linear damping term is included in the equations of motion:



\[

\-\\gamma\\dot{\\theta}\_j

]



where (\\gamma) is the damping coefficient.



With damping disabled, the idealized system conserves mechanical energy up to numerical integration error.



With damping enabled, mechanical energy is dissipated over time.



\---



\## Energy Analysis



The simulation calculates the three main energy quantities:



\### Kinetic Energy



\[

T=

\\frac{1}{2}

\\sum\_{i=1}^{N}

m\_i

\\left(

\\dot{x}\_i^2+\\dot{y}\_i^2

\\right)

]



\### Potential Energy



\[

V=

\\sum\_{i=1}^{N}m\_i g y\_i

]



\### Total Mechanical Energy



\[

E=T+V

]



These quantities are plotted during the simulation to provide a direct way of examining energy conservation and dissipation.



\---



\## Initial Conditions



The simulation provides multiple initial configurations.



\### Horizontal



The pendulum links are initially placed in a horizontal configuration.



\### Inverted



The links are initially placed in an inverted configuration.



\### Random



The angular coordinates are initialized with a random configuration, allowing different nonlinear trajectories to be explored.



\---



\## Features



\### Interactive Pendulum Visualization



The simulation renders the motion of the entire N-link system on a browser canvas.



\### Variable Number of Links



The number of pendulum links can be changed to study how increasing system complexity affects the dynamics.



\### Adjustable Gravity



The gravitational acceleration can be modified to investigate its influence on the system.



\### Damping Control



The damping coefficient can be changed to compare dissipative and approximately conservative motion.



\### Angle Plots



The angular displacement of each link can be plotted as a function of time.



\### Energy Plots



Kinetic, potential, and total mechanical energy can be visualized simultaneously.



\### Animation Controls



The interface includes:



\* Play/pause

\* Reset

\* Animation speed control

\* Zoom and pan

\* Trail visualization

\* Fullscreen charts



\---



\## Web Interface



The project uses \*\*Flask\*\* to provide a lightweight web interface.



The simulation is requested from the backend through the `/api/simulate` endpoint.



A typical simulation request contains parameters such as:



```json

{

&#x20;   "N": 12,

&#x20;   "g": 9.81,

&#x20;   "damping": 0.005,

&#x20;   "preset": "horizontal",

&#x20;   "duration": 15

}

```



The backend performs the numerical integration and returns the simulation data to the frontend.



The frontend then uses the returned data to animate the pendulum and generate the plots.



\---



\## Project Structure



```text

N-Link-Pendulum-Simulation/

│

├── app.py

├── README.md

├── requirements.txt

├── .gitignore

│

└── screenshots/

&#x20;   ├── simulation.png

&#x20;   ├── angles.png

&#x20;   └── energy.png

```



> The screenshot files are optional and can be added after the project is running.



\---



\## Installation



\### 1. Clone the repository



```bash

git clone https://github.com/YOUR-USERNAME/N-Link-Pendulum-Simulation.git

cd N-Link-Pendulum-Simulation

```



\### 2. Create a virtual environment



It is recommended to use a virtual environment:



```bash

python -m venv .venv

```



Activate it on Windows:



```cmd

.venv\\Scripts\\activate

```



\### 3. Install dependencies



```bash

pip install -r requirements.txt

```



\---



\## Running the Simulation



Start the Flask server:



```bash

python app.py

```



The application will run locally at:



```text

http://127.0.0.1:5000

```



Open this address in a web browser to access the simulation.



\---



\## Requirements



The project uses:



\* \*\*Python\*\*

\* \*\*NumPy\*\* — numerical array and matrix operations

\* \*\*SciPy\*\* — numerical integration of the equations of motion

\* \*\*Flask\*\* — backend web application

\* \*\*Chart.js\*\* — interactive plotting

\* \*\*MathJax\*\* — rendering mathematical equations in the interface



Python dependencies are listed in `requirements.txt`.



\---



\## Computational Workflow



The overall simulation pipeline is:



```text

User Parameters

&#x20;     │

&#x20;     ▼

Initial Conditions

&#x20;     │

&#x20;     ▼

N-Link Equations of Motion

&#x20;     │

&#x20;     ▼

Mass Matrix M(θ)

&#x20;     │

&#x20;     ▼

Solve Mα = b

&#x20;     │

&#x20;     ▼

Angular Accelerations

&#x20;     │

&#x20;     ▼

RK45 Numerical Integration

&#x20;     │

&#x20;     ▼

θ(t), θ̇(t)

&#x20;     │

&#x20;     ├───────────────┐

&#x20;     ▼               ▼

Link Positions     Energy

&#x20;     │               │

&#x20;     ▼               ▼

Animation         Graphs

```



\---



\## Results and Observations



The simulation can be used to investigate several characteristics of coupled nonlinear systems.



\### Coupled Motion



The motion of one pendulum link affects the motion of the other links through the configuration-dependent mass matrix and nonlinear coupling terms.



\### Increasing System Complexity



Increasing (N) introduces additional degrees of freedom and stronger opportunities for energy transfer between modes of motion.



\### Energy Behavior



For negligible damping, total mechanical energy should remain approximately constant, with deviations primarily associated with numerical integration error.



For nonzero damping, the total mechanical energy decreases as energy is dissipated.



\### Sensitivity to Initial Conditions



The system can exhibit strongly different trajectories for different initial configurations, particularly as the number of links increases.



\---



\## Screenshots



Add screenshots of the running simulation here.



\### Simulation



\### Angular Displacement



\### Energy



\---



\## Limitations



The current model makes several simplifying assumptions:



\* Pendulum links are treated using prescribed point masses at their ends.

\* Link lengths and masses are configurable but the current frontend uses default values.

\* Air resistance and other complex fluid effects are not modeled.

\* The damping model is simplified as a linear angular damping term.

\* Numerical accuracy depends on the integration tolerances and system configuration.

\* Very large values of (N) increase the computational cost of repeatedly solving the mass-matrix system.



\---



\## Future Improvements



Possible extensions include:



\* \[ ] Separate frontend HTML, CSS, and JavaScript from the Flask backend

\* \[ ] Add independent mass and length controls for every link

\* \[ ] Implement additional numerical integration methods

\* \[ ] Compare RK45 with symplectic integrators for long-term energy behavior

\* \[ ] Add phase-space plots

\* \[ ] Add Poincaré sections

\* \[ ] Quantify sensitivity to initial conditions

\* \[ ] Estimate Lyapunov exponents

\* \[ ] Add parameter-sweep functionality

\* \[ ] Improve performance for large (N)

\* \[ ] Add 3D visualization

\* \[ ] Export simulation data as CSV

\* \[ ] Compare numerical results with analytical solutions for small-angle motion



\---



\## Why This Project?



The N-link pendulum provides a compact example of several important concepts in computational physics and engineering:



\* Classical mechanics

\* Lagrangian mechanics

\* Nonlinear dynamics

\* Coupled differential equations

\* Numerical integration

\* Matrix methods

\* Energy conservation

\* Computational visualization

\* Scientific programming



The project therefore serves as both a physical model and a computational framework for exploring nonlinear dynamical systems.



\---



\## References



The theoretical formulation is based on standard treatments of analytical mechanics and nonlinear dynamics.



Suggested references:



1\. H. Goldstein, C. Poole, and J. Safko, \*Classical Mechanics\*.

2\. L. D. Landau and E. M. Lifshitz, \*Mechanics\*.

3\. S. H. Strogatz, \*Nonlinear Dynamics and Chaos\*.

4\. SciPy documentation — `scipy.integrate.solve\_ivp`.



\---



\## License



This project is intended for educational and research purposes.



You may add a specific open-source license such as the MIT License if you intend to allow others to freely reuse and modify the code.



