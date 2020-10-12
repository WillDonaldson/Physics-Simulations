from numpy import cos, sin, pi, sqrt
from time import time
import matplotlib.pyplot as plt 
import matplotlib.animation as animation

fig = plt.figure()
ax1 = fig.add_subplot(111, aspect='equal')
ax1.set_xlim(-1, 1)
ax1.set_ylim(-1.5, 0.5)
ax1.set_title("Pendulum Animation")
ax1.grid()
line, = ax1.plot([], [], 'o-', lw=1)

theta_0 = pi/8      #the method used in this program is only valid for theta_0 << 1 [rad]
g = 9.81
L = 1.0
omega = sqrt(g/L)
origin = (0,0)
dt = 0.05
t_elapsed = 0.0

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