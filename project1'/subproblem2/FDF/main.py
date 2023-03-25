from utlis import *
import matplotlib.pyplot as plt

a = 10

# Initial condition
y0 = np.array([1, 0]) # molecules/cm^3


# Define the RHS of the ODE
def rhs(t: int , y: list[float]) -> np.array:
    y1, y2 = y
    dy1dt = y2
    dy2dt = -1*a*a * y1
    return np.array([
        dy1dt,
        dy2dt
    ])
    
    
    
# Time span
t_span = (0, 2*np.pi) # seconds
# timestep
k = dt = 0.01 # seconds

########################
f'''Analysis convergence'''
z = a *k

########################


# define the exact solution
exact_y1 = lambda t: np.cos(a*t)
exact_y2 = lambda t: -a * np.sin(a*t)




if __name__ == "__main__":
    
    solver = 'Forward Euler1'
    ode = ODESolver(rhs, dt, t_span, y0, solver=solver)
    ode.solve()

    t_fe = np.arange(t_span[0],t_span[1],dt)

    

    # plot the results for y1
    fig_1, axs_1 = plt.subplots(1, 2,
                                figsize=(10,3),
                                constrained_layout=True)
    axs_1 = axs_1.ravel()
    #plt.figure()
    axs_1[0].plot(np.arange(ode.n_steps+1) * dt, ode.sol[:,0], label=f'{solver}')
    axs_1[0].plot(t_fe, exact_y1(t_fe), label='Exact')
    axs_1[0].set_title(f'y1(t) - BDF1 vs. Exact, k = {k}; |1+z| = {abs(1+z)}')
    axs_1[0].set_xlabel('t')
    axs_1[0].set_ylabel('y1(t)')
    axs_1[0].legend()


    # plot the results for y2
    axs_1[1].plot(np.arange(ode.n_steps+1) * dt, ode.sol[:,1], label=f'{solver}')
    axs_1[1].plot(t_fe, exact_y2(t_fe), label='Exact')
    axs_1[1].set_title(f'y2(t) - BDF1 vs. Exact, k = {k}; |1+z| = {abs(1+z)}')
    axs_1[1].set_xlabel('t')
    axs_1[1].set_ylabel('y2(t)')
    axs_1[1].legend()

    plt.savefig(f'{solver}_k={k}.png')
    plt.show()