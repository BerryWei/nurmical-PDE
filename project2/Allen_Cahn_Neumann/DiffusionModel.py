import numpy as np
from scipy.sparse import coo_matrix
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.sparse
import scipy.sparse.linalg
from etc import read_yaml
from tqdm import tqdm

class DiffusionSimulation:
    def __init__(self, dic):
        self.setting = dic['setting']
        self._initialize_parameters()
        self._create_grid()
        self._build_interior_boundary_grid()
        self._build_laplacian_matrix()
        self.UUs = []
        
        # initial condition
        self.u0 = dic['u0']
        
    def _initialize_parameters(self):
        self.Lx = self.setting['Geometry']['Lx']
        self.Ly = self.setting['Geometry']['Ly']
        self.Tstop = self.setting['Time']['Tstop']
        self.dt = self.setting['Time']['dt']
        self.Nx = self.setting['Mesh']['Nx']
        self.Ny = self.setting['Mesh']['Ny']
        self.alpha = self.setting['diff_coef']
        self.dx = self.Lx / self.Nx
        self.dy = self.Ly / self.Ny
        self.xmin = -self.Lx / 2
        self.ymin = -self.Ly / 2
        self.xmax = self.xmin + self.Lx
        self.ymax = self.ymin + self.Ly

    def _create_grid(self):
        x_array = np.linspace(self.xmin, self.xmax, self.Nx)
        y_array = np.linspace(self.ymin, self.ymax, self.Ny)
        self.XX, self.YY = np.meshgrid(x_array, y_array)
        self.II, self.JJ = np.meshgrid(np.arange(self.Nx), np.arange(self.Ny))
        self.X = self.XX.flatten()
        self.Y = self.YY.flatten()
        self.I = self.II.flatten()
        self.J = self.JJ.flatten()

    def _build_interior_boundary_grid(self):
        isbdy = np.logical_or.reduce([self.I == 0, self.I == self.Nx - 1, self.J == 0, self.J == self.Ny - 1])
        isint = np.logical_not(isbdy)
        self.interior2G = int2G = np.where(isint)[0]
        self.bdy2G = bdy2G = np.where(isbdy)[0]
        self.bdy_left = self.sub2ind([self.Ny, self.Nx], self.I[bdy2G], self.J[bdy2G] - 1)
        self.bdy_down = self.sub2ind([self.Ny, self.Nx], self.I[bdy2G] - 1, self.J[bdy2G])
        self.bdy_right = self.sub2ind([self.Ny, self.Nx], self.I[bdy2G], self.J[bdy2G] + 1)

        self.int_up = self.sub2ind([self.Ny, self.Nx], self.I[int2G] + 1, self.J[int2G])
        self.int_down = self.sub2ind([self.Ny, self.Nx], self.I[int2G] - 1, self.J[int2G])
        self.int_left = self.sub2ind([self.Ny, self.Nx], self.I[int2G], self.J[int2G] - 1)
        self.int_right = self.sub2ind([self.Ny, self.Nx], self.I[int2G], self.J[int2G] + 1)
        self.n_interior = len(self.interior2G)
        self.n_bdy = len(self.bdy2G)
        self.n_grid = self.Nx * self.Ny

    def _build_laplacian_matrix(self):
        # for Diricklet boundary condition --> interior grid
        rows = np.concatenate([np.array(range(self.n_interior)) for _ in range(5)])
        cols = np.hstack([self.interior2G, self.int_left, self.int_right, self.int_up, self.int_down])
        vals = np.hstack([( -2 / self.dx ** 2 - 2 / self.dy ** 2) * np.ones(self.n_interior),
                        1 / self.dx ** 2 * np.ones(self.n_interior),
                        1 / self.dx ** 2 * np.ones(self.n_interior),
                        1 / self.dy ** 2 * np.ones(self.n_interior),
                        1 / self.dy ** 2 * np.ones(self.n_interior)])
        
        self.L_IG = coo_matrix((vals, (rows, cols)), shape=(self.n_interior, self.n_grid)).tocsr()
        self.L_II = self.L_IG.tocsc()[:, self.interior2G]

    def _build_laplacian_matrix_Nuemman(self):
        # for Nuemman boundary condition --> all grid
        f'''TO DO'''
        

    def sub2ind(self, shape, I, J):
        return I * shape[1] + J
    



    def run_simulation_Diricklet(self):
        U = self.u0(self.X, self.Y).astype(float)
        isbdy = np.logical_or.reduce([self.I == 0, self.I == self.Nx - 1, self.J == 0, self.J == self.Ny - 1])
        U[isbdy] = 0
        self.UUs.append(U.copy().reshape([self.Nx, self.Ny]))

        M = scipy.sparse.eye(self.n_interior)
        n_dt = int(np.floor(self.Tstop / self.dt))

        for _ in tqdm(range(n_dt)):
            U[self.interior2G] = scipy.sparse.linalg.spsolve(
                M - self.dt * self.alpha * self.L_II,
                M * U[self.interior2G]
            )
            self.UUs.append(U.copy().reshape([self.Nx, self.Ny]))


    

    def display_animation(self):
        fig, ax = plt.subplots(figsize=(6, 6))
        cax = ax.imshow(self.UUs[0], cmap='viridis', extent=(self.xmin, self.xmax, self.ymin, self.ymax), origin='lower')
        cbar = fig.colorbar(cax)
        ax.set_title("Diffusion Simulation")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_xlim(self.xmin, self.xmax)
        ax.set_ylim(self.ymin, self.ymax)

        time_text = ax.text(0.05, 0.95, '', transform=ax.transAxes, color='white')

        def update(frame):
            cax.set_data(self.UUs[frame])
            time_text.set_text(f't = {frame * self.dt:.2f}')
            return cax, time_text,

        ani = animation.FuncAnimation(fig, update, frames=len(self.UUs), interval=60, blit=True)
        plt.show()

    
    def display_animation(self, save_as='animation.mp4'):
        import matplotlib.animation as animation
        from matplotlib.animation import FFMpegWriter
        fig, ax = plt.subplots(figsize=(6, 6))
        cax = ax.imshow(self.UUs[0], cmap='viridis', extent=(self.xmin, self.xmax, self.ymin, self.ymax), origin='lower')
        cbar = fig.colorbar(cax)
        ax.set_title("Diffusion Simulation")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_xlim(self.xmin, self.xmax)
        ax.set_ylim(self.ymin, self.ymax)

        time_text = ax.text(0.05, 0.95, '', transform=ax.transAxes, color='white')

        def update(frame):
            cax.set_data(self.UUs[frame])
            time_text.set_text(f't = {frame * self.dt:.2f}')
            return cax, time_text,

        ani = animation.FuncAnimation(fig, update, frames=len(self.UUs), interval=60, blit=True)

        if save_as is not None:
            writer = FFMpegWriter(fps=60, bitrate=1800)
            ani.save(save_as, writer=writer)

        plt.show()
        
       


