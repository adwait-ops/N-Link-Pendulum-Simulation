import numpy as np
from flask import Flask, jsonify, request, render_template_string
from scipy.integrate import solve_ivp

app = Flask(__name__)

def derive_n_pendulum_equations(t, y, N, m, L, g, damping):
    theta = y[:N]
    omega = y[N:]
    
    idx = np.maximum.outer(np.arange(N), np.arange(N))
    m_sum = np.sum(m) - np.cumsum(m)[idx] + m[idx]
    
    d_theta = theta[:, None] - theta[None, :]
    M = m_sum * L[:, None] * L[None, :] * np.cos(d_theta)
    
    centrifugal = np.sum(m_sum * L[:, None] * L[None, :] * np.sin(d_theta) * (omega**2)[None, :], axis=1)
    m_gravity = np.cumsum(m[::-1])[::-1]
    gravity = m_gravity * g * L * np.sin(theta)
    
    b = -centrifugal - gravity - damping * omega
    alpha = np.linalg.solve(M, b)
    
    return np.concatenate([omega, alpha])

@app.route('/api/simulate', methods=['POST'])
def simulate():
    data = request.get_json(silent=True) or {}
    N = int(data.get('N', 12))
    g = float(data.get('g', 9.81))
    damping = float(data.get('damping', 0.005))
    duration = float(data.get('duration', 15.0))
    preset = data.get('preset', 'horizontal')
    mass_val = float(data.get('mass', 1.0))
    len_val = float(data.get('length', 0.25))
    fps = 60
    
    m = np.ones(N) * mass_val
    L = np.ones(N) * len_val
    
    if preset == 'horizontal':
        initial_theta = np.full(N, np.pi / 2) + np.random.uniform(-0.01, 0.01, N)
    elif preset == 'inverted':
        initial_theta = np.full(N, np.pi - 0.02) + np.random.uniform(-0.01, 0.01, N)
    else:
        initial_theta = np.random.uniform(-np.pi, np.pi, N)
        
    initial_omega = np.zeros(N)
    y0 = np.concatenate([initial_theta, initial_omega])
    
    t_eval = np.linspace(0, duration, int(duration * fps))
    
    sol = solve_ivp(
        fun=derive_n_pendulum_equations,
        t_span=(0, duration),
        y0=y0,
        t_eval=t_eval,
        args=(N, m, L, g, damping),
        method='RK45',
        rtol=1e-3,
        atol=1e-3
    )
    
    theta_sol = sol.y[:N, :]
    omega_sol = sol.y[N:, :]
    time_steps = len(sol.t)
    
    positions = []
    kinetic_e = []
    potential_e = []
    
    for t_idx in range(time_steps):
        th = theta_sol[:, t_idx]
        om = omega_sol[:, t_idx]
        
        pts = [[0.0, 0.0]]
        x, y = 0.0, 0.0
        
        for i in range(N):
            x += L[i] * np.sin(th[i])
            y -= L[i] * np.cos(th[i])
            pts.append([x, y])
            
        positions.append(pts)
        
        T_val, V_val = 0.0, 0.0
        for i in range(N):
            vx_i = np.sum([L[k] * np.cos(th[k]) * om[k] for k in range(i+1)])
            vy_i = np.sum([L[k] * np.sin(th[k]) * om[k] for k in range(i+1)])
            y_i = pts[i+1][1]
            
            T_val += 0.5 * m[i] * (vx_i**2 + vy_i**2)
            V_val += m[i] * g * y_i
            
        kinetic_e.append(float(T_val))
        potential_e.append(float(V_val))
        
    total_e = (np.array(kinetic_e) + np.array(potential_e)).tolist()
    
    return jsonify({
        'time': sol.t.tolist(),
        'angles': theta_sol.tolist(),
        'positions': positions,
        'length_val': len_val,
        'energy': {
            'kinetic': kinetic_e,
            'potential': potential_e,
            'total': total_e
        },
        'N': N
    })

@app.route('/')
def index():
    return render_template_string(HTML_TEMPLATE)

HTML_TEMPLATE = r"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>N-Link Pendulum Simulation</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <script>
        MathJax = {
            tex: {
                inlineMath: [['$', '$'], ['\\(', '\\)']],
                displayMath: [['$$', '$$'], ['\\[', '\\]']]
            },
            svg: { fontCache: 'global' }
        };
    </script>
    <script type="text/javascript" id="MathJax-script" async
        src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js">
    </script>

    <style>
        * { box-sizing: border-box; margin: 0; padding: 0; }
        body {
            background-color: #0b0b0e;
            color: #d1d5db;
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Arial, sans-serif;
            padding: 24px;
            min-height: 100vh;
            overflow-x: hidden;
        }
        
        .header-title {
            font-size: 1.5rem;
            font-weight: 600;
            color: #f3f4f6;
            margin-bottom: 20px;
            border-bottom: 1px solid #27272a;
            padding-bottom: 10px;
            max-width: 1600px;
            margin-left: auto;
            margin-right: auto;
        }

        .main-layout {
            display: grid;
            grid-template-columns: minmax(0, 1fr) minmax(0, 1fr);
            gap: 24px;
            max-width: 1600px;
            margin: 0 auto;
            align-items: start;
            width: 100%;
        }

        @media (max-width: 1100px) {
            .main-layout { grid-template-columns: minmax(0, 1fr); }
        }

        .panel {
            background: #121215;
            border: 1px solid #27272a;
            border-radius: 8px;
            padding: 16px;
            display: flex;
            flex-direction: column;
            gap: 16px;
            min-width: 0;
            width: 100%;
        }

        .canvas-container {
            position: relative;
            width: 100%;
            height: 480px;
            background-color: #000000;
            border-radius: 6px;
            border: 1px solid #27272a;
            overflow: hidden;
        }

        canvas#simCanvas {
            width: 100%;
            height: 100%;
            display: block;
            cursor: grab;
        }
        canvas#simCanvas:active { cursor: grabbing; }

        .overlay-tools {
            position: absolute;
            top: 10px;
            right: 10px;
            display: flex;
            gap: 6px;
            z-index: 5;
        }

        .btn-sm {
            background: #27272a;
            color: #e4e4e7;
            border: 1px solid #3f3f46;
            padding: 4px 10px;
            border-radius: 4px;
            font-size: 0.8rem;
            cursor: pointer;
            transition: background 0.15s ease;
        }
        .btn-sm:hover { background: #3f3f46; }

        .action-bar {
            display: flex;
            gap: 10px;
        }
        .btn-action {
            flex: 1;
            background: #27272a;
            color: #f4f4f5;
            border: 1px solid #3f3f46;
            padding: 10px;
            border-radius: 6px;
            font-size: 0.9rem;
            font-weight: 500;
            cursor: pointer;
        }
        .btn-action:hover { background: #3f3f46; }
        
        .btn-primary {
            background: #2563eb;
            color: #ffffff;
            border: none;
            padding: 10px;
            border-radius: 6px;
            font-weight: 600;
            cursor: pointer;
            grid-column: 1 / -1;
            transition: background 0.15s ease;
        }
        .btn-primary:hover { background: #1d4ed8; }
        .btn-primary:disabled { background: #1e3a8a; opacity: 0.6; cursor: wait; }

        .controls-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
            gap: 12px;
            background: #09090c;
            padding: 14px;
            border-radius: 6px;
            border: 1px solid #27272a;
        }

        .field-group {
            display: flex;
            flex-direction: column;
            gap: 4px;
            font-size: 0.82rem;
            color: #a1a1aa;
        }

        .field-header {
            display: flex;
            justify-content: space-between;
        }

        input, select {
            background: #18181b;
            color: #f4f4f5;
            border: 1px solid #27272a;
            padding: 6px 10px;
            border-radius: 4px;
            font-size: 0.85rem;
            outline: none;
        }
        input:focus, select:focus { border-color: #3b82f6; }

        .chart-card {
            background: #09090c;
            border: 1px solid #27272a;
            border-radius: 6px;
            padding: 12px;
            display: flex;
            flex-direction: column;
            gap: 8px;
            position: relative;
            box-sizing: border-box;
            width: 100%;
            min-width: 0;
        }

        .chart-card:fullscreen {
            width: 100vw;
            height: 100vh;
            padding: 24px;
            background: #09090c;
            justify-content: center;
        }

        .chart-header {
            display: flex;
            justify-content: space-between;
            align-items: center;
        }

        .chart-title {
            font-size: 0.85rem;
            font-weight: 600;
            color: #d4d4d8;
        }

        .chart-box {
            position: relative;
            height: 210px;
            width: 100%;
            min-width: 0;
        }

        .chart-card:fullscreen .chart-box {
            height: calc(100vh - 90px);
        }

        details {
            background: #09090c;
            border: 1px solid #27272a;
            border-radius: 6px;
            padding: 12px 16px;
            font-size: 0.88rem;
            line-height: 1.65;
            color: #a1a1aa;
            min-width: 0;
        }

        summary {
            font-weight: 600;
            cursor: pointer;
            color: #38bdf8;
            outline: none;
        }

        .explanation-body {
            margin-top: 10px;
            padding-top: 10px;
            border-top: 1px solid #27272a;
        }

        .explanation-body h4 {
            color: #f3f4f6;
            margin: 12px 0 4px 0;
            font-size: 0.92rem;
        }

        .bottom-section {
            max-width: 1600px;
            margin: 24px auto 0 auto;
            width: 100%;
        }

        .math-block {
            background: #030305;
            padding: 12px;
            border-radius: 6px;
            margin: 10px 0;
            overflow-x: auto;
            border: 1px solid #27272a;
            color: #f4f4f5;
        }

        .theory-content h3 {
            color: #ffffff;
            margin: 20px 0 8px 0;
            font-size: 1.05rem;
            border-bottom: 1px solid #27272a;
            padding-bottom: 4px;
        }
    </style>
</head>
<body>

    <div class="header-title">N-Link Pendulum Simulation</div>

    <div class="main-layout">
        <!-- Left Panel: Simulation Canvas & Controls -->
        <div class="panel">
            <div class="canvas-container">
                <div class="overlay-tools">
                    <button class="btn-sm" onclick="zoomCanvas(1.2)" title="Zoom In">+</button>
                    <button class="btn-sm" onclick="zoomCanvas(0.8)" title="Zoom Out">−</button>
                    <button class="btn-sm" onclick="resetCanvasTransform()">Reset View</button>
                </div>
                <canvas id="simCanvas" width="800" height="480"></canvas>
            </div>

            <div class="action-bar">
                <button id="playBtn" class="btn-action" onclick="togglePlay()">Pause</button>
                <button class="btn-action" onclick="resetAnimation()">Reset Trajectory</button>
            </div>

            <div class="controls-grid">
                <div class="field-group">
                    <label>Links (N):</label>
                    <input type="number" id="linkCount" value="12" min="1" max="20">
                </div>
                <div class="field-group">
                    <label>Initial Configuration:</label>
                    <select id="preset">
                        <option value="horizontal">Horizontal Drop</option>
                        <option value="inverted">Inverted Drop</option>
                        <option value="chaos">Random Drop</option>
                    </select>
                </div>
                <div class="field-group">
                    <label>Trail Color:</label>
                    <select id="trailStyle">
                        <option value="rainbow">Rainbow</option>
                        <option value="cyan">Cyan</option>
                        <option value="fire">Thermal</option>
                        <option value="matrix">Green</option>
                        <option value="monochrome">White</option>
                    </select>
                </div>
                <div class="field-group">
                    <div class="field-header">
                        <label>Speed ($\Delta t$):</label>
                        <span id="speedVal">1.00x</span>
                    </div>
                    <input type="range" id="simSpeed" min="0.05" max="2.0" step="0.05" value="1.0" oninput="updateSpeedLabel(this.value)">
                </div>
                <div class="field-group">
                    <label>Damping ($\gamma$):</label>
                    <input type="number" id="damping" value="0.005" step="0.002">
                </div>
                <div class="field-group">
                    <label>Gravity ($g$):</label>
                    <input type="number" id="gravity" value="9.81" step="0.5">
                </div>
                <div class="field-group">
                    <label>Trail Length:</label>
                    <input type="range" id="trailLen" min="50" max="1000" value="450">
                </div>
                <button id="runBtn" class="btn-primary" onclick="runSimulation()">Re-run Simulation</button>
            </div>
        </div>

        <!-- Right Panel: Dynamic Graphs & Explanations -->
        <div class="panel">
            <div style="display: flex; justify-content: space-between; align-items: center;">
                <span style="font-weight: 600; color: #f4f4f5; font-size: 0.95rem;">System</span>
                <div style="display: flex; gap: 4px;">
                    <button class="btn-sm" onclick="setChartWindow(3)">3s</button>
                    <button class="btn-sm" onclick="setChartWindow(6)">6s</button>
                    <button class="btn-sm" onclick="setChartWindow(15)">15s</button>
                </div>
            </div>

            <!-- Angle Graph Card -->
            <div class="chart-card" id="angleCard">
                <div class="chart-header">
                    <span class="chart-title">Link Angles $\theta_1 \dots \theta_N$ (rad)</span>
                    <button class="btn-sm" id="angleBtn" onclick="toggleFullscreen('angleCard', 'angleBtn')">Fullscreen</button>
                </div>
                <div class="chart-box"><canvas id="angleChart"></canvas></div>
            </div>

            <details>
                <summary>Angle Graph Details</summary>
                <div class="explanation-body">
                    <h4>Coordinate Reference Frame</h4>
                    <p>The graph plots the absolute spatial inclination angle $\theta_k(t)$ of each link $k \in \{1, \dots, N\}$ in radians, relative to the downward vertical gravitational vector ($0\text{ rad} = 6\text{ o'clock}$).</p>
                    
                    <h4>Key Behavior Indicators</h4>
                    <ul style="margin-left: 18px; margin-top: 6px;">
                        <li><b>Harmonic Oscillations ($\pm\pi$ boundary):</b> Smooth wave trajectories staying within $\approx -3.14$ and $+3.14\text{ rad}$ indicate links swinging back and forth without completing a full loop.</li>
                        <li><b>Rotational Inversion & Continuous Accumulation:</b> When a line climbs or drops linearly across multiples of $2\pi$, that link has acquired enough kinetic energy to continuously spin $360^\circ$ around its joint.</li>
                        <li><b>Phase Divergence & Mode Coupling:</b> Higher-indexed bobs ($\theta_N$) exhibit high-frequency jitter superimposed on low-frequency waves, illustrating how energy cascades from top links down to the tip.</li>
                    </ul>
                </div>
            </details>

            <!-- Energy Graph Card -->
            <div class="chart-card" id="energyCard">
                <div class="chart-header">
                    <span class="chart-title">Mechanical Energy Conservation (J)</span>
                    <button class="btn-sm" id="energyBtn" onclick="toggleFullscreen('energyCard', 'energyBtn')">Fullscreen</button>
                </div>
                <div class="chart-box"><canvas id="energyChart"></canvas></div>
            </div>

            <details>
                <summary>Energy Graph Details</summary>
                <div class="explanation-body">
                    <h4>Energy Partition Breakdown</h4>
                    <ul style="margin-left: 18px; margin-top: 6px;">
                        <li><b style="color: #f87171;">Kinetic Energy ($T$):</b> Scalar sum of translational velocities across all $N$ masses: $T = \frac{1}{2}\sum m_i (\dot{x}_i^2 + \dot{y}_i^2)$.</li>
                        <li><b style="color: #4ade80;">Potential Energy ($V$):</b> Gravitational potential defined relative to origin $(0,0)$: $V = g \sum m_i y_i$.</li>
                        <li><b style="color: #facc15;">Total Energy ($E = T + V$):</b> Global Hamiltonian scalar tracking energy conservation.</li>
                    </ul>

                    <h4>Dissipation Mechanics</h4>
                    <p style="margin-top: 6px;">When Rayleigh damping $\gamma = 0$, Total Energy $E(t)$ forms a strictly horizontal invariant line. When $\gamma > 0$, $E(t)$ exhibits monotonic non-reversible exponential decay due to viscous torque friction $-\gamma \dot{\theta}_j$.</p>
                </div>
            </details>
        </div>
    </div>

    <!-- Bottom Section: Full Mathematical Physics Theory -->
    <div class="bottom-section">
        <details>
            <summary>Mathematics & Non-Linear Dynamics Theory</summary>
            <div class="explanation-body theory-content">
                <h3>1. Generalized Kinematics</h3>
                <p>An $N$-link planar pendulum constrained to a 2D vertical plane has $N$ degrees of freedom completely described by the angular vector $\boldsymbol{\theta}(t) = [\theta_1(t), \theta_2(t), \dots, \theta_N(t)]^T$, where $\theta_k$ measures the inclination of link $k$ relative to the downward vertical axis. Position vectors $\mathbf{r}_i = (x_i, y_i)$ of each bob $i \in \{1, \dots, N\}$ are given recursively from the fixed origin $(0,0)$:</p>
                <div class="math-block">
                    $$x_i(\boldsymbol{\theta}) = \sum_{k=1}^i L_k \sin\theta_k, \quad y_i(\boldsymbol{\theta}) = -\sum_{k=1}^i L_k \cos\theta_k$$
                </div>
                <p>Differentiating with respect to time yields the velocity vectors $\mathbf{v}_i = \dot{\mathbf{r}}_i = (\dot{x}_i, \dot{y}_i)$:</p>
                <div class="math-block">
                    $$\dot{x}_i = \sum_{k=1}^i L_k \dot{\theta}_k \cos\theta_k, \quad \dot{y}_i = \sum_{k=1}^i L_k \dot{\theta}_k \sin\theta_k$$
                </div>

                <h3>2. Hamiltonian Formulation & Canonical Phase Space</h3>
                <p>The total scalar kinetic energy $T$ and potential energy $V$ are expressed as:</p>
                <div class="math-block">
                    $$T(\boldsymbol{\theta}, \dot{\boldsymbol{\theta}}) = \frac{1}{2}\sum_{i=1}^N m_i (\dot{x}_i^2 + \dot{y}_i^2) = \frac{1}{2} \dot{\boldsymbol{\theta}}^T M(\boldsymbol{\theta}) \dot{\boldsymbol{\theta}}$$
                </div>
                <div class="math-block">
                    $$V(\boldsymbol{\theta}) = g \sum_{i=1}^N m_i y_i = -g \sum_{i=1}^N \left( \sum_{k=i}^N m_k \right) L_i \cos\theta_i$$
                </div>

                <h3>3. Exact Matrix Formulation of Equations of Motion</h3>
                <p>Applying the Euler-Lagrange equations with Rayleigh dissipation function $R = \frac{1}{2}\gamma \sum_{j=1}^N \dot{\theta}_j^2$ yields the matrix relation $M(\boldsymbol{\theta})\ddot{\boldsymbol{\theta}} = \mathbf{b}(\boldsymbol{\theta}, \dot{\boldsymbol{\theta}})$:</p>
                <div class="math-block">
                    $$M_{j,k}(\boldsymbol{\theta}) = \left(\sum_{i=\max(j,k)}^{N} m_i\right) L_j L_k \cos(\theta_j - \theta_k)$$
                </div>
                <div class="math-block">
                    $$b_j(\boldsymbol{\theta}, \dot{\boldsymbol{\theta}}) = -\sum_{k=1}^N \left(\sum_{i=\max(j,k)}^N m_i\right) L_j L_k \sin(\theta_j - \theta_k) \dot{\theta}_k^2 - \left(\sum_{i=j}^N m_i\right) g L_j \sin\theta_j - \gamma \dot{\theta}_j$$
                </div>
            </div>
        </details>
    </div>

    <script>
        const canvas = document.getElementById('simCanvas');
        const ctx = canvas.getContext('2d');
        
        let angleChart = null;
        let energyChart = null;
        let animationFrameId = null;
        let isPlaying = true;
        
        let simData = null;
        let currentFrameFloat = 0;
        let tipHistory = [];
        
        // Canvas Interaction State
        let zoomScale = 1.0;
        let panX = 0;
        let panY = 0;
        let isDragging = false;
        let dragDistance = 0;
        let startX = 0, startY = 0;

        let chartWindowSec = 15;

        function updateSpeedLabel(val) {
            document.getElementById('speedVal').innerText = parseFloat(val).toFixed(2) + 'x';
        }

        // Canvas Navigation Listeners
        canvas.addEventListener('wheel', (e) => {
            e.preventDefault();
            const factor = e.deltaY < 0 ? 1.15 : 0.85;
            zoomCanvas(factor);
        });

        canvas.addEventListener('mousedown', (e) => {
            isDragging = true;
            dragDistance = 0;
            startX = e.clientX - panX;
            startY = e.clientY - panY;
        });

        window.addEventListener('mousemove', (e) => {
            if (!isDragging) return;
            dragDistance += Math.hypot(e.movementX, e.movementY);
            panX = e.clientX - startX;
            panY = e.clientY - startY;
        });

        window.addEventListener('mouseup', () => {
            if (isDragging && dragDistance < 5) {
                togglePlay();
            }
            isDragging = false;
        });

        function zoomCanvas(factor) {
            zoomScale = Math.min(Math.max(0.2, zoomScale * factor), 8.0);
            if (!isPlaying) drawFrame();
        }

        function resetCanvasTransform() {
            zoomScale = 1.0;
            panX = 0;
            panY = 0;
            if (!isPlaying) drawFrame();
        }

        function toggleFullscreen(cardId, btnId) {
            const card = document.getElementById(cardId);
            const btn = document.getElementById(btnId);
            
            if (!document.fullscreenElement) {
                card.requestFullscreen().then(() => {
                    btn.innerText = 'Exit Fullscreen';
                }).catch(err => alert(err.message));
            } else {
                document.exitFullscreen().then(() => {
                    btn.innerText = 'Fullscreen';
                });
            }
        }

        document.addEventListener('fullscreenchange', () => {
            if (!document.fullscreenElement) {
                document.getElementById('angleBtn').innerText = 'Fullscreen';
                document.getElementById('energyBtn').innerText = 'Fullscreen';
            }
            
            const chartCanvases = document.querySelectorAll('.chart-box canvas');
            chartCanvases.forEach(c => {
                c.style.width = '100%';
                c.style.height = '100%';
            });

            setTimeout(() => {
                if (angleChart) angleChart.resize();
                if (energyChart) energyChart.resize();
                drawFrame();
            }, 100);
        });

        function setChartWindow(sec) {
            chartWindowSec = sec;
            if (!isPlaying) drawFrame();
        }

        async function runSimulation() {
            const runBtn = document.getElementById('runBtn');
            runBtn.innerText = 'Calculating ODE...';
            runBtn.disabled = true;

            try {
                const N = parseInt(document.getElementById('linkCount').value) || 12;
                const damping = parseFloat(document.getElementById('damping').value) || 0.005;
                const g = parseFloat(document.getElementById('gravity').value) || 9.81;
                const preset = document.getElementById('preset').value;

                const res = await fetch('/api/simulate', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ N, damping, g, preset, duration: 15 })
                });
                
                const newData = await res.json();
                
                // Safely update charts BEFORE assigning new global simData
                initCharts(newData);
                simData = newData;
                
                currentFrameFloat = 0;
                tipHistory = [];
                
                if (!animationFrameId) animate();
            } catch (e) {
                console.error("Simulation error:", e);
                alert("Failed to compute simulation. Check terminal for server logs.");
            } finally {
                runBtn.innerText = 'Simulate';
                runBtn.disabled = false;
            }
        }

        function togglePlay() {
            isPlaying = !isPlaying;
            document.getElementById('playBtn').innerText = isPlaying ? 'Pause' : 'Play';
        }

        function resetAnimation() {
            currentFrameFloat = 0;
            tipHistory = [];
            drawFrame();
        }

        function animate() {
            if (simData && simData.positions) {
                if (isPlaying) {
                    const speedMultiplier = parseFloat(document.getElementById('simSpeed').value) || 1.0;
                    currentFrameFloat += speedMultiplier;
                    
                    if (currentFrameFloat >= simData.positions.length) {
                        currentFrameFloat = 0;
                        tipHistory = [];
                    }
                }
                drawFrame();
            }
            animationFrameId = requestAnimationFrame(animate);
        }

        function drawFrame() {
            if (!simData || !simData.positions) return;

            let frameIdx = Math.floor(currentFrameFloat);
            if (frameIdx >= simData.positions.length) frameIdx = 0;

            ctx.clearRect(0, 0, canvas.width, canvas.height);
            
            const lengthVal = simData.length_val || 0.25;
            const totalReach = simData.N * lengthVal;
            const baseScale = totalReach > 0 ? (canvas.height * 0.38) / totalReach : 100;
            const scale = baseScale * zoomScale;
            
            const cx = canvas.width / 2 + panX;
            const cy = canvas.height * 0.35 + panY;

            const framePos = simData.positions[frameIdx];
            const tip = framePos[framePos.length - 1];
            
            const tipX = cx + tip[0] * scale;
            const tipY = cy - tip[1] * scale;
            
            const maxHist = parseInt(document.getElementById('trailLen').value) || 450;
            
            if (isPlaying) {
                tipHistory.push({x: tipX, y: tipY});
                if (tipHistory.length > maxHist) tipHistory.shift();
            }

            const trailStyle = document.getElementById('trailStyle').value;

            // Render Trail
            for (let i = 1; i < tipHistory.length; i++) {
                ctx.beginPath();
                ctx.moveTo(tipHistory[i-1].x, tipHistory[i-1].y);
                ctx.lineTo(tipHistory[i].x, tipHistory[i].y);
                
                const alpha = i / tipHistory.length;
                
                if (trailStyle === 'rainbow') {
                    const hue = (i * 1.8 + frameIdx) % 360;
                    ctx.strokeStyle = `hsla(${hue}, 100%, 60%, ${alpha})`;
                } else if (trailStyle === 'cyan') {
                    ctx.strokeStyle = `rgba(0, 240, 255, ${alpha})`;
                } else if (trailStyle === 'fire') {
                    const hue = (i * 0.45);
                    ctx.strokeStyle = `hsla(${hue}, 100%, 50%, ${alpha})`;
                } else if (trailStyle === 'matrix') {
                    ctx.strokeStyle = `rgba(34, 197, 94, ${alpha})`;
                } else {
                    ctx.strokeStyle = `rgba(255, 255, 255, ${alpha * 0.9})`;
                }
                
                ctx.lineWidth = 2 * Math.sqrt(zoomScale);
                ctx.stroke();
            }

            // Render Pendulum Rods
            ctx.strokeStyle = 'rgba(220, 225, 240, 0.8)';
            ctx.lineWidth = 1.8 * Math.sqrt(zoomScale);
            ctx.beginPath();
            for (let i = 0; i < framePos.length; i++) {
                const px = cx + framePos[i][0] * scale;
                const py = cy - framePos[i][1] * scale;
                if (i === 0) ctx.moveTo(px, py);
                else ctx.lineTo(px, py);
            }
            ctx.stroke();

            // Render Bobs
            for (let i = 0; i < framePos.length; i++) {
                const px = cx + framePos[i][0] * scale;
                const py = cy - framePos[i][1] * scale;
                ctx.beginPath();
                
                if (i === 0) {
                    ctx.arc(px, py, 4 * Math.sqrt(zoomScale), 0, 2 * Math.PI);
                    ctx.fillStyle = '#64748b';
                } else if (i === framePos.length - 1) {
                    ctx.arc(px, py, 6 * Math.sqrt(zoomScale), 0, 2 * Math.PI);
                    ctx.fillStyle = '#ffffff';
                } else {
                    ctx.arc(px, py, 3 * Math.sqrt(zoomScale), 0, 2 * Math.PI);
                    ctx.fillStyle = '#a1a1aa';
                }
                ctx.fill();
            }

            updateCharts(frameIdx);
        }

        function initCharts(data) {
            const ctxA = document.getElementById('angleChart').getContext('2d');
            const ctxE = document.getElementById('energyChart').getContext('2d');
            
            if (angleChart) angleChart.destroy();
            if (energyChart) energyChart.destroy();

            const angleDatasets = [];
            const colors = ['#f472b6', '#a78bfa', '#38bdf8', '#4ade80', '#facc15', '#fb923c'];

            for (let i = 0; i < data.N; i++) {
                angleDatasets.push({
                    label: `θ${i+1}`,
                    data: [],
                    borderColor: colors[i % colors.length],
                    borderWidth: 1.2,
                    pointRadius: 0
                });
            }

            const chartOptions = {
                responsive: true,
                maintainAspectRatio: false,
                animation: false,
                plugins: { legend: { display: false } },
                scales: {
                    x: {
                        type: 'linear',
                        ticks: { color: '#71717a', callback: (v) => v.toFixed(1) + 's' },
                        grid: { color: '#27272a' }
                    },
                    y: { ticks: { color: '#71717a' }, grid: { color: '#27272a' } }
                }
            };

            angleChart = new Chart(ctxA, {
                type: 'line',
                data: { datasets: angleDatasets },
                options: JSON.parse(JSON.stringify(chartOptions))
            });

            const energyOptions = JSON.parse(JSON.stringify(chartOptions));
            energyOptions.plugins.legend = { display: true, labels: { color: '#a1a1aa', boxWidth: 10, font: { size: 10 } } };

            energyChart = new Chart(ctxE, {
                type: 'line',
                data: {
                    datasets: [
                        { label: 'Kinetic', data: [], borderColor: '#f87171', borderWidth: 1.2, pointRadius: 0 },
                        { label: 'Potential', data: [], borderColor: '#4ade80', borderWidth: 1.2, pointRadius: 0 },
                        { label: 'Total', data: [], borderColor: '#facc15', borderWidth: 1.5, pointRadius: 0 }
                    ]
                },
                options: energyOptions
            });
        }

        function updateCharts(frameIdx) {
            if (!angleChart || !energyChart || !simData || !simData.angles) return;

            const currentTime = simData.time[frameIdx];
            const endIdx = frameIdx + 1;

            let minX = 0;
            let maxX = 15.0;
            
            if (chartWindowSec < 15) {
                minX = Math.max(0, currentTime - chartWindowSec);
                maxX = minX + chartWindowSec;
            }

            // Bound iteration by both dataset length and current angles array length
            const datasetCount = Math.min(angleChart.data.datasets.length, simData.angles.length);
            for (let i = 0; i < datasetCount; i++) {
                if (!simData.angles[i]) continue;
                const pts = [];
                for (let k = 0; k < endIdx; k++) {
                    pts.push({ x: simData.time[k], y: simData.angles[i][k] });
                }
                angleChart.data.datasets[i].data = pts;
            }

            angleChart.options.scales.x.min = minX;
            angleChart.options.scales.x.max = maxX;
            angleChart.update('none');

            const kPts = [], pPts = [], tPts = [];
            for (let k = 0; k < endIdx; k++) {
                const t = simData.time[k];
                kPts.push({ x: t, y: simData.energy.kinetic[k] });
                pPts.push({ x: t, y: simData.energy.potential[k] });
                tPts.push({ x: t, y: simData.energy.total[k] });
            }

            energyChart.data.datasets[0].data = kPts;
            energyChart.data.datasets[1].data = pPts;
            energyChart.data.datasets[2].data = tPts;

            energyChart.options.scales.x.min = minX;
            energyChart.options.scales.x.max = maxX;
            energyChart.update('none');
        }

        window.onload = runSimulation;
    </script>
</body>
</html>"""

if __name__ == '__main__':
    app.run(debug=True, port=5000)
