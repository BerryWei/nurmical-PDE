from utlis import *
import matplotlib.pyplot as plt


# define reaction rates
f_k1 = lambda t: 1e-2 * max(0, np.sin(2 * np.pi * t / (24*60*60))) # s^-1, t_d = 24*60*60 s
k2 = 1e5 # s^-1
k3 = 1e-16 # cm^3 molecule^-1 s^-1

# Initial condition
y0 = np.array([0, 0, 5e11, 8e11]) # molecules/cm^3


# Define the RHS of the ODE
def rhs(t: int , y: list[float]) -> np.array:
    c1, c2, c3, c4 = y
    k1 = f_k1(t) # s^-1, t_d = 24*60*60 s
    return np.array([
        k1 * c3 - k2 * c1,
        k1 * c3 - k3 * c2 * c4,
        k3 * c2 * c4 - k1 * c3,
        k2 * c1 - k3 * c2 * c4
    ])
    
# Time span
t_span = (0, 2*24*60*60) # seconds
# timestep
dt = 60 # seconds




if __name__ == "__main__":
    
    ode = ODESolver(rhs, dt, t_span, y0, solver='newton')
    ode.solve()

    #plot
    

    plt.plot(np.arange(ode.n_steps+1) * dt / (24*60*60), ode.sol[:,0], label='O')
    plt.plot(np.arange(ode.n_steps+1) * dt / (24*60*60), ode.sol[:,1], label='NO')
    plt.plot(np.arange(ode.n_steps+1) * dt / (24*60*60), ode.sol[:,2], label='NO2')
    plt.plot(np.arange(ode.n_steps+1) * dt / (24*60*60), ode.sol[:,3], label='O3')
    plt.legend()
    plt.xlabel('Time (days)')
    plt.ylabel('Concentration (molecules/cm^3)')
    plt.show()
    plt.savefig('solution.png')