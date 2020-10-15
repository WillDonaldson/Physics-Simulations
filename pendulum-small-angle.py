from numpy import cos, sin, pi
import numpy as np
from time import time
import matplotlib.pyplot as plt 
import matplotlib.animation as animation
from scipy.integrate import solve_ivp

fig = plt.figure()
ax1 = fig.add_subplot(121, aspect='equal', xlim = (-1, 1), ylim = (-1.5, 0.5), title = "Pendulum Animation")
ax2 = fig.add_subplot(122, xlim = (-2*pi, 2*pi), ylim = (-10, 10), title = "Phase Space Plot")
ax2.set_xlabel(r"$\Theta$[rad]")
ax2.set_ylabel(r"$\dot{\Theta}$[rad/s]")

ax1.grid()
ax2.grid()

line, = ax1.plot([], [], 'o-', lw=1)
point, = ax2.plot([],[], 'ro')

theta_0 = [pi/8, 0]    # theta_0[1] = initial angular velocity
g = 9.81
L = 1.0
m = 1.0
omega = np.sqrt(g/L)
origin = (0,0)
dt = 0.05
t_elapsed = 0.0
t_span = [0,30]

def Hamiltonian(q, p):
    #Calculates the Hamiltonian of a simple pendulum
    H = p**2 / (2*m*L**2) + m*g*L*(1-cos(q))
    return H

#Create points for plotting the phase space
x = np.arange(-2*pi, 2*pi, 0.05)
y = np.arange(-10, 10, 0.05)
Theta, Theta_dot = np.meshgrid(x, y)            # poor notation, needs renaming. This variable is different to the theta in calculate_position()
q = Theta
p = m * L**2 * Theta_dot
cs = ax2.contour(Theta, Theta_dot, Hamiltonian(q,p))
ts = np.linspace(t_span[0], t_span[1], 600)

def eqn(t, theta_0):
    # f = [theta, theta_dot]
    # returns f'
    return [theta_0[1], -omega**2 * sin(theta_0[0])]

pendulum_state = solve_ivp(eqn, t_span, theta_0, t_eval = ts)



def calculate_position(i):
    #theta = theta_0 * cos(omega * t)            # valid for small angles, theta_0 << 1 [rad]
    #theta = odeint(f, theta_0, t, args = omega)
    #theta = solve_ivp(eqn, t_span, theta_0, t_eval = ts)
    #theta_dot = -theta_0[0] * omega * sin(omega * t) + theta_0[1]
    theta = pendulum_state.y[0][i]
    theta_dot = pendulum_state.y[1][i] 
    xData = [origin[0], L * sin(theta)]
    yData = [origin[1], -L * cos(theta)]
    return theta, theta_dot, xData, yData

def animate(i):
    global t_elapsed
    t_elapsed += dt
    #theta, theta_dot, x, y = calculate_position(t_elapsed)
    theta, theta_dot, x, y = calculate_position(i)
    line.set_data(x, y)   
    point.set_data(theta, theta_dot)
    return line, point,

t0 = time()
animate(1)                          #sample time to evaluate function
t1 = time()
interval = 1000 * dt - (t1 - t0)

ani = animation.FuncAnimation(fig, animate, frames = 100, interval = interval)
plt.show()