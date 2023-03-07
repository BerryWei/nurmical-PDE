import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

a = -100

# Initial condition
y0 = np.array([1]) # molecules/cm^3


# Define the RHS of the ODE
def rhs(t: int , y: list[float]) -> np.array:
    return np.array([
        a * y[0]
    ])
    
# Time span
t_span = (0, 1) # seconds
# timestep
k = dt = 5e-1 # seconds
t_eval = np.arange(t_span[0], t_span[1], dt) # create the time vector for evaluation

########################
f'''Analysis convergence'''
z = a *k
########################

# Define the exact solution
def y_exact(t,y0):
    return y0*np.exp(a*t)

solver_names = ['RK45', 'RK23', 'DOP853', 'Radau', 'BDF', 'LSODA']


if __name__ == "__main__":


    t_fe = np.arange(t_span[0],t_span[1],1e-4)
    y_exact_vals = y_exact(t_fe, y0)

    #plot
    plt.figure()
    for i in range(len(solver_names)):
        solver_name = solver_names[i]
        solver = solve_ivp(rhs, t_span, y0, method=solver_name, t_eval=t_eval)
        plt.plot(solver.t, solver.y[0], label=solver_name) # plot the solution for the current solver

    plt.plot(t_fe, y_exact_vals, 'k-', label='Exact')
    plt.legend()
    plt.xlabel('t')
    plt.ylabel('y(t)')
    plt.title(f'a = {a}, dt = {dt:.1e}, z = {z}')
    plt.savefig(f'solution1_dt={dt:.1e}.png')
    plt.show()






    

    plt.plot(np.arange(ode.n_steps+1) * dt, ode.sol[:,0], label=solver)
    plt.plot(t_fe, y_exact_vals, 'k-', label='Exact')
    plt.legend()
    plt.xlabel('t')
    plt.ylabel('y(t)')
    plt.title(f'a = {a}, k = {k}, |1+z| = {abs(1+z)}')
    
    plt.savefig(f'solution1_k={k:1.0e}.png')
    plt.show()