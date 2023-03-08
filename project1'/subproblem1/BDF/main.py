from utlis import *
import matplotlib.pyplot as plt

a = -100

# Initial condition
y0 = np.array([1]) # molecules/cm^3


# Define the RHS of the ODE
def rhs(t: int , y: list[float]) -> np.array:
    return np.array([
        a * y
    ])
    
# Time span
t_span = (0, 1) # seconds
# timestep
k = dt = 1e-1 # seconds

########################
f'''Analysis convergence'''
z = a *k

########################

# Define the exact solution
def y_exact(t,y0):
    return y0*np.exp(a*t)





if __name__ == "__main__":
    
    solver = 'BDF1-newton'
    ode = ODESolver(rhs, dt, t_span, y0, solver=solver)
    ode.solve()

    t_fe = np.arange(t_span[0],t_span[1],dt)
    y_exact_vals = y_exact(t_fe, y0)

    #plot

    plt.plot(np.arange(ode.n_steps+1) * dt, ode.sol[:,0], label=solver)
    plt.plot(t_fe, y_exact_vals, 'k-', label='Exact')
    plt.legend()
    plt.xlabel('t')
    plt.ylabel('y(t)')
    plt.title(f'a = {a}, k = {k}, |1-z| = {abs(1-z)}')
    
    plt.savefig(f'solution1_k={k:1.0e}.png')
    plt.show()