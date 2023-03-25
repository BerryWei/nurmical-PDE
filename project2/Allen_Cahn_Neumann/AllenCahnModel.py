from Model import PDESimulation
from scipy.sparse import coo_matrix
import scipy
import numpy as np
from tqdm import tqdm


class AllenCahnSimulation(PDESimulation):
    def __init__(self, dic):
        super().__init__(dic)
        self._build_laplacian_matrix_neumann()
        

    def _initialize_parameters(self):
        super()._initialize_parameters()
        # self.dt = 0.1 * self.dx ** 2
        self.epsilon = self.setting['epsilon']

        # FOR FFT EXAMPLE
        self.Lx = self.setting['Geometry']['Lx']
        self.Ly = self.setting['Geometry']['Ly']
        self.Tstop = self.setting['Time']['Tstop']
        self.dt = self.setting['Time']['dt']
        self.Nx = self.setting['Mesh']['Nx']
        self.Ny = self.setting['Mesh']['Ny']
        self.dx = self.Lx / self.Nx
        self.dy = self.Ly / self.Ny
        self.xmin = 0
        self.ymin = 0
        self.xmax = self.xmin + self.Lx
        self.ymax = self.ymin + self.Ly

    def run_simulation_IE(self):
        f"""Here, we assumed that flux = 0
        Diffusion term: implicit scheme
        Reaction term: explicit scheme	
        """
        # initial conditon
        U = self.u0(self.X, self.Y).astype(float)

        f'''U[isbdy]: u(x,y,t) on boundary grid'''
        isbdy = np.logical_or.reduce([self.I == 0, self.I == self.Nx - 1, self.J == 0, self.J == self.Ny - 1])
        isint = ~isbdy
        self.UUs.append(U.copy().reshape([self.Nx, self.Ny]))

        I = scipy.sparse.identity(self.Nx * self.Ny)
        I_bdy = scipy.sparse.identity(len(np.where(isbdy == True)[0]))
        I_int = scipy.sparse.identity(len(np.where(isint == True)[0]))
        n_dt = int(np.floor(self.Tstop / self.dt))

        pbar = tqdm(range(n_dt))
        for _ in pbar:
            # diffusion term (Time implicit, space center difference)
            pbar.set_description(f"Processing...")

            U = U.flatten()
            # interior grids
            LHS = I_int - self.dt * self.epsilon * self.L_GG[:, isint][isint, :]
            RHS = U[isint] + self.dt * (U[isint] - U[isint] ** 3) 
            U_int = scipy.sparse.linalg.spsolve(LHS, RHS)
            
            # boundary grids
            LHS = I_bdy - 2 * self.dt * self.epsilon * self.L_GG[:, isbdy][isbdy, :]
            RHS = U[isbdy] + 2 * self.dt * (U[isbdy] - U[isbdy] ** 3) 
            U_bdy = scipy.sparse.linalg.spsolve(LHS, RHS)
            
            # fill the value into U
            U[isint] = U_int
            U[isbdy] = U_bdy

            self.UUs.append(U.copy().reshape([self.Nx, self.Ny]))

    def run_simulation_IE_splitting(self):
        f"""Here, we assumed that flux = 0
        Diffusion term: Implicit scheme
        Reaction term: Explicit scheme
        Solve with splitting method
        """
        # initial conditon
        U = self.u0(self.X, self.Y).astype(float)

        f'''U[isbdy]: u(x,y,t) on boundary grid'''
        isbdy = np.logical_or.reduce([self.I == 0, self.I == self.Nx - 1, self.J == 0, self.J == self.Ny - 1])
        isint = ~isbdy
        self.UUs.append(U.copy().reshape([self.Nx, self.Ny]))

        I = scipy.sparse.identity(self.Nx * self.Ny)
        I_bdy = scipy.sparse.identity(len(np.where(isbdy == True)[0]))
        I_int = scipy.sparse.identity(len(np.where(isint == True)[0]))

        n_dt = int(np.floor(self.Tstop / self.dt))

        pbar = tqdm(range(n_dt))
        for _ in pbar:
            # diffusion term (Time implicit, space center difference)
            pbar.set_description(f"Processing...")

            U = U.flatten()
            term1_bdy = scipy.sparse.linalg.spsolve(I_bdy - 2* self.dt * self.epsilon * self.L_GG[:, isbdy][isbdy, :], U[isbdy] )
            term1_int = scipy.sparse.linalg.spsolve(I_int -    self.dt * self.epsilon * self.L_GG[:, isint][isint, :], U[isint] )
            term1 = np.zeros(self.n_grid)
            term1[isbdy] = term1_bdy
            term1[isint] = term1_int
            # Neumann boundary *2 on boundary grid

            term2_bdy = U[isbdy] + 2*self.dt * (-U[isbdy]**3 + U[isbdy])
            term2_int = U[isint] +   self.dt * (-U[isint]**3 + U[isint])
            term2 = np.zeros(self.n_grid)
            term2[isbdy] = term2_bdy
            term2[isint] = term2_int

            U = term1*term2/U

            self.UUs.append(U.copy().reshape([self.Nx, self.Ny]))

    def run_simulation_IE_fft(self):
        f"""Here, we assumed that flux = 0
        Diffusion term: Implicit scheme
        Reaction term: Explicit scheme
        Solve with FFT
        """
        # initial conditon
        U = self.u0(self.X, self.Y).astype(float).reshape([self.Nx, self.Ny])

        I = complex(0, 1)
        k_x = np.array([I * n for n in list(range(0, self.Nx // 2)) + [0] + list(range(-self.Nx // 2 + 1, 0))])
        k_y = k_x

        kxx = np.zeros((self.Nx, self.Nx), dtype=complex)
        kyy = np.zeros((self.Nx, self.Nx), dtype=complex)
        for i in range(self.Nx):
            for j in range(self.Nx):
                kxx[i, j] = k_x[i] ** 2
                kyy[i, j] = k_y[j] ** 2

        self.UUs.append(U.copy())
        n_dt = int(np.floor(self.Tstop / self.dt))


        pbar = tqdm(range(n_dt))
        for _ in pbar:
            # diffusion term (Time implicit, space center difference)
            pbar.set_description(f"Processing...")

            v_hat = np.fft.fft2(U)            
            # calculate nonlinear term in real space
            v_nl = np.array(U**3 , dtype=complex)
            # FFT for nonlinear and linear terms
            v_nl = np.fft.fft2(v_nl)
            v_hat = (v_hat * ( 1 + 1 / self.dt) - v_nl)
            v_hat = v_hat / (1 / self.dt - (kxx + kyy) * self.epsilon)  # Implicit/Explicit timestepping
            
            U = np.real(np.fft.ifft2(v_hat))

            self.UUs.append(U.copy())

