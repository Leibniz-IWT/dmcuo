import numpy as np

#% Settings
def load_settings():
    #% Options
    class Opts():
        def __init__(self):
           self.t_stop  = 250
           self.errorbars = 1  # % Display errorbars (0 - off = 1 - on)
           #self.N_nodes = 2500  # % Number of nodes in numerical grid for solution of the diffusion equation (diffusion_controlled_release()) in -
           self.N_nodes = 500  # % Number of nodes in numerical grid for solution of the diffusion equation (diffusion_controlled_release()) in -

           self.t_max = 1e-6  # % Max. timestep for solution of diffusion equation in s
           self.disp_intermediate_results = 1  # % Display intermediate results during itteration (0 - off = 1 - on)
           self.disp_particle_profile = 1  # % Display concentration profile in particles (c(r =t)) from diffusion equation (0 - off = 1 - on)

    opts = Opts()
    #% Note: N_nodes  = 2500 required for sufficient accuracy
    
    #% Physical constants and initial parameters
    class Params():
        def __init__(self):
           self.rho_CuO = 6.31 * 1e-21  # % Density of CuO in g/nm�
           self.M_CuO = 79.545 / 1000  # % Molar mass of CuO in g/mmol
           self.N_A = 6.022e23  # % Avogadro constant in #/mol
           self.n_data = 4  # % Number of datasets for evaluation in -
           self.c0 = [2.5, 0, 5, 0]  # % Initial molar concentrations of [CuO = Cu2+ = Amino Acid = Crystals] in mmol
           self.d0 = [11.8, 12.0, 10.3, 10.7]  # % Initial particle diameters as obtained from BET measurements in nm
           self.m_r = 3.5  # % Reaction order in -
           self.n_r = 4  # % Reaction order in -
           self.ratioFe_Cu = [0, 0.01, 0.06, 0.1]  # % Iron-copper ratio in -

    params = Params()
    #% Color definitions
    black = np.array([0, 0, 0])/255.0
    cyan =np.array( [0, 147, 221])/255.0
    green = np.array([132, 194, 37])/255.0
    red = np.array([218, 37, 29])/255.0
    colors = [black, cyan, green, red]

    return opts, params, colors