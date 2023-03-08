import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# define the system of ODEs
def system(t, y, a):
    y1, y2 = y
    dy1dt = y2
    dy2dt = -1*a*a * y1
    return [dy1dt, dy2dt]

# set the initial conditions and value of a
y0 = [1, 0]
a = 10

# define the time span and time step size
t_span = [0, 2*np.pi]
t_step = 0.01

# define the exact solution
exact_y1 = lambda t: np.cos(a*t)
exact_y2 = lambda t: -a * np.sin(a*t)

method = 'LSODA'
# define the RK45 solver
rk45_sol = solve_ivp(lambda t, y: system(t, y, a), t_span, y0, method=method, t_eval=np.arange(t_span[0], t_span[1], t_step))

# plot the results for y1
fig_1, axs_1 = plt.subplots(1, 2,
                            figsize=(10,3),
                            constrained_layout=True)
axs_1 = axs_1.ravel()
#plt.figure()
axs_1[0].plot(rk45_sol.t, rk45_sol.y[0], label='RK45')
axs_1[0].plot(rk45_sol.t, exact_y1(rk45_sol.t), label='Exact')
axs_1[0].set_title(f'Comparison of y1(t) - {method} vs. Exact')
axs_1[0].set_xlabel('t')
axs_1[0].set_ylabel('y1(t)')
axs_1[0].legend()


# plot the results for y2

axs_1[1].plot(rk45_sol.t, rk45_sol.y[1], label='RK45')
axs_1[1].plot(rk45_sol.t, exact_y2(rk45_sol.t), label='Exact')
axs_1[1].set_title(f'Comparison of y2(t) - {method} vs. Exact')
axs_1[1].set_xlabel('t')
axs_1[1].set_ylabel('y2(t)')
axs_1[1].legend()

plt.savefig(f'{method}.png')
plt.show()