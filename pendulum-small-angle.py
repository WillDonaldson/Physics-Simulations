from numpy import cos, sin, pi
import numpy as np
from time import time
import matplotlib.pyplot as plt 
import matplotlib.animation as animation

fig = plt.figure()
ax1 = fig.add_subplot(121, aspect='equal', xlim = (-1, 1), ylim = (-1.5, 0.5), title = "Pendulum Animation")
ax2 = fig.add_subplot(122, xlim = (-2*pi, 2*pi), ylim = (-10, 10), title = "Phase Space Plot")
ax2.set_xlabel(r"$\Theta$[rad]")
ax2.set_ylabel(r"$\dot{\Theta}$[rad/s]")

ax1.grid()
ax2.grid()

line, = ax1.plot([], [], 'o-', lw=1)

theta_0 = pi/8      #the method used in this program is only valid for theta_0 << 1 [rad]
g = 9.81
L = 1.0
m = 1.0
omega = np.sqrt(g/L)
origin = (0,0)
dt = 0.05
t_elapsed = 0.0

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

def calculate_position(t):
    theta = theta_0 * cos(omega * t)            # valid for small angles
    xData = [origin[0], L * sin(theta)]
    yData = [origin[1], -L * cos(theta)]
    return xData, yData

def animate(i):
    global t_elapsed
    t_elapsed += dt
    line.set_data(calculate_position(t_elapsed))   
    return line, 

t0 = time()
print(t_elapsed)
animate(0)                          #sample time to evaluate function
t1 = time()
interval = 1000 * dt - (t1 - t0)

ani = animation.FuncAnimation(fig, animate, frames = 100, interval = interval)
plt.show()