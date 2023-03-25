from AllenCahnModel import AllenCahnSimulation
from etc import read_yaml
import numpy as np
from pathlib import Path
import math

np.random.seed(16)


def IC1(x, y):
    return np.random.rand(len(x))*2-1

def IC2(x, y):
    return np.exp(-((x - 1) ** 2 + (y - 2) ** 2) / 10)


def IC3(x, y):
    return np.sin(2*math.pi*x) + 0.5*np.cos(4*math.pi*y)


if __name__ == '__main__':
    configPath = Path(f'./AllenCahn.yaml')
    config={}
    config['setting'] = read_yaml(configPath)  # Assuming the 'read_yaml' function has been defined
    config['u0'] = IC1



    sim = AllenCahnSimulation(config)
    sim.run_simulation_IE_splitting()
    sim.display_animation_save(
        title = 'Allen Cahn 2D model with Neumann BC(IE_splitting)',
        save_as= 'IE_splitting.mp4'
    )
    #########################################################################
    sim = AllenCahnSimulation(config)
    sim.run_simulation_IE()
    sim.display_animation_save(
        title = 'Allen Cahn 2D model with Neumann BC(IE_coupling)',
        save_as= 'IE_coupling.mp4'
    )
    #########################################################################
    sim = AllenCahnSimulation(config)
    sim.run_simulation_IE_fft()
    sim.display_animation_save(
        title = 'Allen Cahn 2D model with Neumann BC(IE_fft)',
        save_as= 'IE_fft.mp4'
    )
