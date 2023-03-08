import numpy as np
from scipy.optimize import root
import sys 

class ODESolver:
    def __init__(self, rhs, dt, t_span, y0, solver='hybird', jac = None) -> None:
        self.rhs = rhs # rhs of ODE system
        self.solver = solver
        self.dt = dt # time step
        self.t_span = t_span # time span
        self.y0 = y0  # initial condition
        self.jac = jac  # Jacobian function
        # ---- 
        self.n_steps = int((t_span[1] - t_span[0]) / dt)
        self.sol = np.zeros((self.n_steps + 1, len(y0)))
        self.sol[0] = y0

        # dictionary of available solvers
        self.dir_solvers = {
            'BDF1-hybird': root,
            'BDF1-newton': self.newton_raphson,
            'BDF1-moore_penrose': self.moore_penrose,
            'Forward Euler1': self.forward_euler1,
        }

    def solve(self, **kwargs):
        solver_fn = self.dir_solvers[self.solver]
        for i in range(self.n_steps):
            y_prev = self.sol[i]
            t_prev = i * self.dt
            t_next = (i + 1) * self.dt
            # find root
            if self.solver == 'BDF1-newton' or self.solver == 'BDF1-moore_penrose':
                res = solver_fn(self.rhs, t_next, y_prev, self.dt)
                self.sol[i + 1] = res

            elif self.solver == 'Forward Euler1':
                res = solver_fn(self.rhs, t_prev, y_prev, self.dt)
                self.sol[i + 1] = res

            else:

                res = solver_fn(lambda x: x - y_prev - self.dt * self.rhs(t_next, x),
                                y_prev, method='hybr', tol=1e-6, jac=self.jac,
                                options={'maxfev': 100})
                self.sol[i + 1] = res.x

    @staticmethod
    def forward_euler1(fun : callable, t_prev:float, x0 : np.ndarray, k:float ):
        return x0 + k * fun(t_prev, x0)

    @staticmethod
    def newton_raphson(fun : callable, t_next:float, x0 : np.ndarray, k:float, eps=1e-6, max_iter=100):
        r'''Find a root of a vector function.
        Parameters
        ----------
        fun: f(t, y)
        x0: initial guess
        '''
        n = len(x0)
        x = x0.copy()
        F = lambda x: x - x0 - k*fun(t_next, x)

        for i in range(max_iter):
            # compute Jacobian
            J = ODESolver.jacobian(F, x)
            # update x
            #dx = - np.linalg.inv(J) @ fun(x) --> slower
            dx = np.linalg.solve(J, -F(x))
            x = x + dx

            # check if convergence
            if np.linalg.norm(dx) < eps:
                return x
        return RuntimeError(f"newton_raphson failed to converge after {max_iter} iterations.")
    

    def moore_penrose(fun, x0, eps=1e-6, max_iter=100):
        n = len(x0)
        x = x0.copy()

        for i in range(max_iter):
            # compute Jacobian
            J = ODESolver.jacobian(fun, x)
            print(J.T @ J)
            # check if Jacobian is singular
            if np.linalg.cond(J) > 1 / sys.float_info.epsilon:
                # compute the pseudo inverse of J
                J_pinv = np.linalg.inv(J.T @ J) @ J.T
                dx = J_pinv @ -fun(x)
            else:
                dx = np.linalg.solve(J, -fun(x))

            x = x + dx

            # check if convergence
            if np.linalg.norm(dx) < eps:
                return x

        raise RuntimeError(f"newton_raphson failed to converge after {max_iter} iterations.")
    
    @staticmethod
    def jacobian(fun : callable, x : np.ndarray, h=1e-6):
        n = len(x)
        J = np.zeros((n, n))
        x0 = x.copy()
        for i in range(n):
            
            # calculate partial derivatives
            x0[i] += h
            f1 = fun(x0)
            x0[i] -= 2 * h
            f2 = fun(x0)
            x0[i] += h
            df = (f1 - f2) / (2 * h)

            # fill in Jacobian matrix
            J[:, i] = df
        return J
