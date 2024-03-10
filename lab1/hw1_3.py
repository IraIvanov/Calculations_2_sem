import numpy as np
import matplotlib.pyplot as mpl
import math
from scipy.optimize import fsolve

td = 24*60*60
T = 2*td
k2 = 10**(5)
k3 = 10**(-16)

def k1(t):

    return max(0, math.sin(2*math.pi*t/td))/100

def F(c, t):

    C = np.array([k1(t)*c[2] - k2*c[0], k1(t)*c[2] - k3*c[1]*c[3], -k1(t)*c[2] + k3*c[1]*c[3], k2*c[0] - k3*c[1]*c[3]])
    return C

def Eu(c_n, delta_t, t):

    c_n_1 = c_n + delta_t*F(c_n, t)
    return c_n_1

C_0 = np.array([0, 0, 5*(10**11), 8*(10**11)])

def eq(c_n, t, delta_t):

    k_1 = k1(t + delta_t)
    const = tuple(F(c_n, t))
    c_n = tuple(c_n)

    def eq1(p):

        x, y, z, a = p
        x_1 = delta_t*(k_1*z - k2*x + const[0])/2  +c_n[0]
        y_1 = delta_t*(k_1*z - k3*y*a + const[1])/2 + c_n[1]

        return (x_1, y_1, delta_t*(k3*y*a - k_1*z + const[2])/2 + c_n[2], delta_t*(k2*x - k3*y*a + const[3])/2 + c_n[3])

    return eq1 

def K_N(c_0, n):

    delta_t = T/n
    solution = [c_0]
    c_1 = [c_0[0]]
    c_2 = [c_0[1]]
    c_3 = [c_0[2]]
    c_4 = [c_0[3]]
    time = [0]

    for i in range(n + 1):

        c_i = solution[i]
        t_i = time[i]
        t_next = t_i + delta_t
        c_next = Eu(c_i, delta_t, t_i)#fsolve(eq(c_i, t_i, delta_t), c_i)
        solution.append(c_next)
        time.append(t_next)
        c_1.append(c_next[0])
        c_2.append(c_next[1])
        c_3.append(c_next[2])
        c_4.append(c_next[3])

    return solution, c_1, c_2, c_3, c_4, time

n = 10000

sol, c_1, c_2, c_3, c_4, time = K_N(C_0, n)

fig, ax = mpl.subplots()

ax.plot(time, c_1, 'r-', label=('C1'))
mpl.grid()
mpl.xlabel('Время')
mpl.ylabel('Значение в точке')
mpl.title('График решений')
mpl.legend(loc='best', bbox_to_anchor = (1, -0.1))
mpl.tight_layout()
mpl.savefig("images/Euler-c1.png")

fig, ax = mpl.subplots()

ax.plot(time, c_2, 'b-', label=('C2'))
mpl.grid()
mpl.xlabel('Время')
mpl.ylabel('Значение в точке')
mpl.title('График решений')
mpl.legend(loc='best', bbox_to_anchor = (1, -0.1))
mpl.tight_layout()
mpl.savefig("images/Euler-c2.png")

fig, ax = mpl.subplots()

ax.plot(time, c_3, 'g-', label=('C3'))
mpl.grid()
mpl.xlabel('Время')
mpl.ylabel('Значение в точке')
mpl.title('График решений')
mpl.legend(loc='best', bbox_to_anchor = (1, -0.1))
mpl.tight_layout()
mpl.savefig("images/Euler-c3.png")

fig, ax = mpl.subplots()

ax.plot(time, c_4, '-', label=('C4'))
mpl.grid()
mpl.xlabel('Время')
mpl.ylabel('Значение в точке')
mpl.title('График решений')
mpl.legend(loc='best', bbox_to_anchor = (1, -0.1))
mpl.tight_layout()
mpl.savefig("images/Euler-c4.png") 

