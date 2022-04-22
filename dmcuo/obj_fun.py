import numpy as np
import scipy.integrate
from dmcuo.reaction_kinetic_ODE_system import reaction_kinetic_ODE_system, dcdt_r
from dmcuo.diffusion_controlled_release import diffusion_controlled_release

#%% Objective function that calculates the root mean square error (RMSE) between experimental data and model for a set of parameters (params_x[#-1]#%% The model data is obtained from the superimposed solution of the kinetically driven relase (solution 1) and the diffusion limited release (solution 2)
def obj_fun(params_x, data, params, opts, colors):
    # init contianers
    k = np.zeros(len(params_x))
    c = np.zeros(8)  # Correct?
    # results storage contianers:
    sim = []
    sim2 = []
    model = {'t': [],
             'c_cu2p_1': [],
             'c_cu2p_2': [],
             }

    minSquareError = 0.0  #%% Error initialization
    for l in range(params.n_data):
        # #%% Prepare experimental data - Concentration profiles
        data_time = data['time'][l]
        data_c_cu2p = data['mean'][l]
        data_time = np.trim_zeros(data_time, 'b')  # trim trailing zeros from data set
        data_c_cu2p = np.trim_zeros(data_c_cu2p, 'b')  # trim trailing zeros from data set

        # #%% Initial conditions and parameters for the ODE system describing the kinetics (ode45 solver)
        # % Initial particle diameter as obtained from BET measurements
        d0 = params.d0[l]  # % in nm
        # % Initial particle radius
        r0 = 0.5*d0  # % in nm
        
        # #%% Fit parameters for the ODE system
        # % Dissolution rate constants
        k[1-1] = params_x[4-1]*1e-36
        k[2-1] = params_x[5-1]*1e-36
        k[3-1] = params_x[6-1]*1e-36
        k[4-1] = params_x[7-1]*1e-36
        # % Diffusion coefficient
        k[5-1] = params_x[2-1]
        # % Atomic surface concentration per unit surface area in #/nm2
        k[6-1] = params_x[1-1]*1000*params.M_CuO/params.rho_CuO/params.N_A/0.073
        
        # #%% Initial concentrations
        c[1-1] = params_x[3-1]  # % C_CuO (variable) in mM
        c[2-1] = params.c0[2-1]  # % C_Cu2+ (fixed) in mM
        c[3-1] = params.c0[3-1]  # % C_AA (amino acid, fixed) in mM
        
        # #%% Calculated values based on input
        V0 = c[1-1]*params.M_CuO/params.rho_CuO  # % Particle volume in nm2
        n_particle = V0/(4/3*np.pi*r0*r0*r0)  # % Number of particles in -
        a0 = 3*V0/r0  # % Total particle surface area in nm2
        c[4-1] = r0  # % r_CuO - Particle radius in nm
        c[5-1] = a0  # % a_CuO - Particle surface area in nm2
        c[6-1] = V0  # % V_CuO - Particle volume in nm2
        
        # #%% Initial iron-copper and copper-iron ratio
        c[7-1] = params.ratioFe_Cu[l]*c[5-1]/k[6-1]# % f_fe
        c[8-1] = (1-params.ratioFe_Cu[l])*c[5-1]/k[6-1]# % f_Cu
        
        # #%% Set initial conditions for the ODE system
        ic = [c[1-1], c[2-1], c[3-1], c[4-1], c[5-1], c[6-1], c[7-1], c[8-1]]

        # #%% Options for integrator:
        ode_options = {'rtol': 2e-14,
                       'atol': 1e-12 }

        # reaction_kinetic_ODE_system(c,k,params,n_particle,l))
        # NOTE: dcdt_r is a wrapper for reaction_kinetic_ODE_system
        y0 = ic
        sol1 = scipy.integrate.solve_ivp(dcdt_r, (data_time[0], data_time[-1]), y0, method='RK45',
                                         t_eval=data_time,
                                         args=(k, params, n_particle, l),
                                         **ode_options)

        # #%% SOLUTION 2: Numerical solution of the diffusion model defined in diffusion_controlled_release
        r_final = sol1.y[4-1,-1]
        sol2 = diffusion_controlled_release(r_final, k, params, opts, c, l, t_eval=data_time)
        sol2.y = np.nan_to_num(sol2.y)
        if l == 1 - 1:
            sol2.y = np.zeros(len(data_time))  # Should be zero since it is pure CuO dissolves without the diffusion
                                               # controlled release, but force in case
            sol2.t = np.zeros(len(data_time))
        else:
            sol2.y = sol2.y[-1, :]

        # #%% Evaluation of solutions 1 and 2 at the experimental data (times)
        sim.append(sol1.y)
        sim2.append(sol2.y)

        # % Plot intermediate results for each iteration step
        #NOTE: plt.pause(0.01) can do this if needed, however, I would not recommend
        #      this in any simulation loop for performance considerations. - Stefan Endres
       # #%% Save model results for plot
        model['t'].append(sol1.t)
        model['c_cu2p_1'].append(sol1.y)
        sol2y = sol2.y
        model['c_cu2p_2'].append(sol2y)

       # #%% Calculation of the objective function, RMSE  # minSquareError = minSquareError + sum((sim{l}(2,:) + sim2{l}(end,:) - conc_eval{l}).^2)/length(conc_eval{l})
       #NOTE: conc_eval became data_c_cu2p after the coversion
        try:
            minSquareError = np.linalg.norm(sol1.y + sol2y - data_c_cu2p) / len(data_c_cu2p)
        except ValueError:
            minSquareError = np.inf
    return [minSquareError, model, sol2]
