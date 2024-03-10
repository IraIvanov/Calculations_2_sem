import numpy as np
import matplotlib.pyplot as mpl
import math

T = 7.416298709205487
y_0 = np.array([1, 0])

def func(y):

    f_y = [0, 0]
    f_y[0] = y[1]
    f_y[1] = -(y[0])*(y[0])*(y[0])
    f_y = np.array(f_y)

    return f_y

def R_K(x_0, n):

    delta_t = T/n
    solution = [x_0]
    u = [x_0[0]]
    v = [x_0[1]]
    time = [0]

    for i in range(n + 1):

        x_i = solution[i]
        t_i = time[i]
        X0 = x_i
        T_0 = t_i
        X1 = X0 + delta_t/2*func(X0)
        T_1 = T_0 + delta_t/2
        X2 = X1 + delta_t/2*func(X1)
        T_2 = T_1 + delta_t/2
        X3 = X2 + delta_t*func(X2)
        T_3 = T_2 + delta_t

        x_next = x_i + delta_t*(func(X0)/6 + func(X1)/3 + func(X2)/3 + func(X3)/6)
        t_next = t_i + delta_t
        u.append(x_next[0])
        v.append(x_next[1])
        solution.append(x_next)
        time.append(t_next)

    return solution, u, v, time

n = 1000
sol, u, v, time = R_K(y_0, n)

mpl.plot(time, u, 'r-', label=('Рещение x'))
mpl.plot(time, v, '-', label=('Рещение y'))
mpl.grid()
mpl.xlabel('Время')
mpl.ylabel('Значение в точке')
mpl.title('График решений')
mpl.legend(loc='best', bbox_to_anchor = (1, -0.1))
mpl.tight_layout()
mpl.show()
#mpl.savefig('images/solution1.png')

n = 1000
norms = []
ns = []

while n < 100000:
    sol, u, v, time = R_K(y_0, n)
    y_T = sol[-1]
    Y = y_T - y_0
    norm = math.sqrt(Y[1]*Y[1] + Y[0]*Y[0])
    ns.append(T/n)
    norms.append(norm)
    n *= 10

mpl.plot(ns, norms, 'ro-', label=('погрешность'))
mpl.grid()
mpl.xlabel('Шаг')
mpl.ylabel('Норма погешности')
mpl.title('График погрешности')
mpl.legend(loc='best', bbox_to_anchor = (1, -0.1))
mpl.tight_layout()
mpl.show()

k, b = np.polyfit(ns, norms, 1)
print(k)

z = 1 +1j

def region(z):

    return abs(z*z*z*z/24 + z*z*z/6 + z*z/2 + z + 1)

x_l = []
y_l = []

x = -5
while x < 10:
    y = -5
    while y < 5:
        z = complex(x,y)
        if region(z) < 1:
            x_l.append(x)
            y_l.append(y)
        y += 0.01
    x += 0.01

mpl.plot(x_l, y_l, '-', label=("Область устойчивости"))
mpl.grid()
mpl.xlabel('Re')
mpl.ylabel('Im')
mpl.legend(loc='best', bbox_to_anchor = (1, -0.1))
mpl.tight_layout()
mpl.savefig('images/region')
