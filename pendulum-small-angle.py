# notes to self
# come back and add in option for rod with non-zero mass
# cleanup seperate function for phase space data
# repair timing issue introduced by t_span

import numpy as np
from time import time
import matplotlib.pyplot as plt 
import matplotlib.animation as animation
from scipy.integrate import solve_ivp

fig = plt.figure()
ax1 = fig.add_subplot(121, aspect='equal', xlim = (-1, 1), ylim = (-1.5, 0.5), title = "Pendulum Animation")
ax2 = fig.add_subplot(122, xlim = (-2*np.pi, 2*np.pi), ylim = (-10, 10), title = "Phase Space Plot")
ax2.set_xlabel(r"$\Theta$[rad]")
ax2.set_ylabel(r"$\dot{\Theta}$[rad/s]")

ax1.grid()
ax2.grid()

line, = ax1.plot([], [], 'o-', lw=1)    # pendulum arm 
point, = ax2.plot([],[], 'ro')          # position in phase space

theta_0 = [np.pi/8, 0]                  # theta_0[1] = initial angular velocity
g = 9.81
L = 1.0
m = 1.0
omega = np.sqrt(g/L)
origin = (0,0)
dt = 0.05
t_span = [0,30]

def Hamiltonian(q, p):
    # calculates the Hamiltonian of a simple pendulum
    H = p**2 / (2*m*L**2) + m*g*L*(1-np.cos(q))
    return H

# create points for plotting
x = np.arange(-2*np.pi, 2*np.pi, 0.05)
y = np.arange(-10, 10, 0.05)
Theta, Theta_dot = np.meshgrid(x, y)            # poor notation, needs renaming. This variable is different to the theta in calculate_position()
q = Theta
p = m * L**2 * Theta_dot
cs = ax2.contour(Theta, Theta_dot, Hamiltonian(q,p))
ts = np.linspace(t_span[0], t_span[1], 600)

def eqn(t, theta_0):
    # f = [theta, theta_dot]
    # returns f'
    return [theta_0[1], -omega**2 * np.sin(theta_0[0])]

pendulum_state = solve_ivp(eqn, t_span, theta_0, t_eval = ts)

def calculate_position(i):
    #theta = theta_0 * cos(np.omega * t)            # valid for small angles, theta_0 << 1 [rad]
    #theta_dot = -theta_0[0] * omega * np.sin(omega * t) + theta_0[1]
    theta = pendulum_state.y[0][i]
    theta_dot = pendulum_state.y[1][i] 
    xData = [origin[0], L * np.sin(theta)]
    yData = [origin[1], -L * np.cos(theta)]
    return theta, theta_dot, xData, yData

def animate(i):
    theta, theta_dot, x, y = calculate_position(i)
    line.set_data(x, y)   
    point.set_data(theta, theta_dot)
    return line, point,

t0 = time()
animate(0)                          #sample time to evaluate function
t1 = time()
interval = 1000 * dt - (t1 - t0)

ani = animation.FuncAnimation(fig, animate, frames = 100, interval = interval)
plt.show()