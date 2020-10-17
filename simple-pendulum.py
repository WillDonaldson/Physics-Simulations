# This program creates an animation of a rod-pendulum (uniform mass) 

# Future improvements could include:
#  - User input for initial conditions
#  - Options for the type of pendulum: simple, rod, compound, and/or N-body pendulum
#  - Dampening factor in the equations of motion

import numpy as np
from time import time
import matplotlib.pyplot as plt 
import matplotlib.animation as animation
from scipy.integrate import solve_ivp

# setup plots
fig = plt.figure()
ax1 = fig.add_subplot(121, aspect='equal', xlim = (-1, 1), ylim = (-1.5, 0.5), title = "Pendulum Animation")
ax2 = fig.add_subplot(122, xlim = (-2*np.pi, 2*np.pi), ylim = (-15, 15), title = "Phase Space Plot")
ax2.set_xlabel(r"$\Theta$[rad]")
ax2.set_ylabel(r"$\dot{\Theta}$[rad/s]")
ax1.grid()
ax2.grid()
line, = ax1.plot([], [], 'o-', lw=1)    # pendulum arm 
point, = ax2.plot([],[], 'ro')          # position in phase space

# pendulum parameters
theta_0 = [np.pi/8, 0.0]                # theta_0[1] = initial angular velocity
g = 9.81
L = 1.0
m = 1.0
I = m*L**2/3                            # moment of inertia for a rod pendulum
omega = np.sqrt((m*g*L)/(2*I))

# animation parameters
origin = [0.0, 0.0]
dt = 0.05
frames = 600
t_span = [0.0, frames * dt]

def Hamiltonian(q, p):
    H = p**2 / (6*m*L**2) + m*g*L/2*(1-np.cos(q))
    return H

def eqn(t, theta_0):
    # f = [theta, theta_dot]
    # returns f'
    return [theta_0[1], -omega**2 * np.sin(theta_0[0])]

ts = np.linspace(t_span[0], t_span[1], frames)
pendulum_state = solve_ivp(eqn, t_span, theta_0, t_eval = ts)

# phase space data points
# (optional) this code snippet could be refactored in terms of pendulum_state.y[][]
x = np.linspace(-2*np.pi, 2*np.pi, frames)
y = np.linspace(-15, 15, frames)
ThetaGrid, Theta_dotGrid = np.meshgrid(x, y)            
q = ThetaGrid                           # generalised coordinate
p = m * L**2 * Theta_dotGrid            # generalise momementum 
cs = ax2.contour(ThetaGrid, Theta_dotGrid, Hamiltonian(q,p))

def animate(i):
    theta = pendulum_state.y[0][i]
    theta_dot = pendulum_state.y[1][i] 
    x = [origin[0], L * np.sin(theta)]
    y = [origin[1], -L * np.cos(theta)]
    line.set_data(x, y)   
    point.set_data(theta, theta_dot)
    return line, point,

t0 = time()
animate(0)                          #sample time required to evaluate function
t1 = time()
interval = 1000 * dt - (t1 - t0)

ani = animation.FuncAnimation(fig, animate, frames = frames, interval = interval)
plt.show()