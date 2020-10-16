# notes to self
# come back and add in option for rod with non-zero mass
# cleanup seperate function for phase space data
# user input for IC, and type of pendulum: simple, rod, compound, add in dampening?
# double pendulum case
# above features are more physics than programming practice, so I may omit them for now

import numpy as np
from time import time
import matplotlib.pyplot as plt 
import matplotlib.animation as animation
from scipy.integrate import solve_ivp

# setup plots
fig = plt.figure()
ax1 = fig.add_subplot(121, aspect='equal', xlim = (-1, 1), ylim = (-1.5, 0.5), title = "Pendulum Animation")
ax2 = fig.add_subplot(122, xlim = (-2*np.pi, 2*np.pi), ylim = (-10, 10), title = "Phase Space Plot")
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
omega = np.sqrt(g/L)

# animation parameters
origin = [0.0, 0.0]
dt = 0.05
frames = 600
t_span = [0.0, frames * dt]

def Hamiltonian(q, p):
    # calculates the Hamiltonian of a simple pendulum
    H = p**2 / (2*m*L**2) + m*g*L*(1-np.cos(q))
    return H

def eqn(t, theta_0):
    # f = [theta, theta_dot]
    # returns f'
    return [theta_0[1], -omega**2 * np.sin(theta_0[0])]

ts = np.linspace(t_span[0], t_span[1], frames)
pendulum_state = solve_ivp(eqn, t_span, theta_0, t_eval = ts)

# phase space data points
# this code snippet could be refactored in terms of pendulum_state.y[][]
# elected to keep using np.meshgrid() for ease of readability
x = np.linspace(-2*np.pi, 2*np.pi, frames)
y = np.linspace(-10, 10, frames)
ThetaGrid, Theta_dotGrid = np.meshgrid(x, y)            
q = ThetaGrid
p = m * L**2 * Theta_dotGrid            # conjugate momementum 
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