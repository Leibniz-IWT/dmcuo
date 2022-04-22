import numpy as np
import scipy.integrate
from dmcuo.reaction_kinetic_ODE_system import reaction_kinetic_ODE_system, dcdt_r
from dmcuo.diffusion_controlled_release import diffusion_controlled_release


# %% Objective function that calculates the root mean square error (RMSE) between experimental data and model for a set of parameters (params_x[#-1]#%% The model data is obtained from the superimposed solution of the kinetically driven relase (solution 1) and the diffusion limited release (solution 2)
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
    # #%% Error initialization
    minSquareError = 0.0
    print(f'data = {data}')
    print(f'opts = {opts}')
    # #%% Model calculation for pure CuO, #1%Fe+CuO, #6%Fe+CuO and 1#0%Fe+CuO
    for l in range(params.n_data):
        # #%% Prepare experimental data - Concentration profiles
        # % Indicies of the data (table with varying length)
        #    first_idx = 1
        #    last_idx = find(data.time(l,:) ~= 0,1,'last')
        # % Remove empty datapoints (zeros) in the table
        #   time_eval{l} = data.time(l,first_idx:last_idx)
        #   conc_eval{l} = data.c_cu2p(l,first_idx:last_idx)
        # % Remove datapoints exceeding opts.t_stop limit  #TODO:
        #  time_eval{l} = time_eval{l}(find(time_eval{l} <= opts.t_stop))
        #   conc_eval{l} = conc_eval{l}(find(time_eval{l} <= opts.t_stop))

        data_time = data['time'][l]
        data_c_cu2p = data['mean'][l]
        data_time = np.trim_zeros(data_time, 'b')  # trim trailing zeros from data set
        data_c_cu2p = np.trim_zeros(data_c_cu2p, 'b')  # trim trailing zeros from data set

        # #%% Initial conditions and parameters for the ODE system describing the kinetics (ode45 solver)
        # % Initial particle diameter as obtained from BET measurements
        d0 = params.d0[l]  # % in nm
        # % Initial particle radius
        r0 = 0.5 * d0  # % in nm

        # #%% Fit parameters for the ODE system
        # % Dissolution rate constants
        k[1 - 1] = params_x[4 - 1] * 1e-36
        k[2 - 1] = params_x[5 - 1] * 1e-36
        k[3 - 1] = params_x[6 - 1] * 1e-36
        k[4 - 1] = params_x[7 - 1] * 1e-36
        # % Diffusion coefficient
        k[5 - 1] = params_x[2 - 1]
        # % Atomic surface concentration per unit surface area in #/nm2
        k[6 - 1] = params_x[1 - 1] * 1000 * params.M_CuO / params.rho_CuO / params.N_A / 0.073

        # #%% Initial concentrations
        c[1 - 1] = params_x[3 - 1]  # % C_CuO (variable) in mM
        c[2 - 1] = params.c0[2 - 1]  # % C_Cu2+ (fixed) in mM
        c[3 - 1] = params.c0[3 - 1]  # % C_AA (amino acid, fixed) in mM

        # #%% Calculated values based on input
        V0 = c[1 - 1] * params.M_CuO / params.rho_CuO  # % Particle volume in nm2
        n_particle = V0 / (4 / 3 * np.pi * r0 * r0 * r0)  # % Number of particles in -
        a0 = 3 * V0 / r0  # % Total particle surface area in nm2
        c[4 - 1] = r0  # % r_CuO - Particle radius in nm
        c[5 - 1] = a0  # % a_CuO - Particle surface area in nm2
        c[6 - 1] = V0  # % V_CuO - Particle volume in nm2

        # #%% Initial iron-copper and copper-iron ratio
        c[7 - 1] = params.ratioFe_Cu[l] * c[5 - 1] / k[6 - 1]  # % f_fe
        c[8 - 1] = (1 - params.ratioFe_Cu[l]) * c[5 - 1] / k[6 - 1]  # % f_Cu

        # #%% Set initial conditions for the ODE system
        ic = [c[1 - 1], c[2 - 1], c[3 - 1], c[4 - 1], c[5 - 1], c[6 - 1], c[7 - 1], c[8 - 1]]

        # #%% Options for ode45
        # ode_opts = odeset('RelTol',2e-14,'AbsTol',1e-12,'Nonnegative',[])
        ode_options = {'rtol': 2e-14,
                       'atol': 1e-12}
        # #%% Timespan of model simulation
        # tSpan = [0:1:opts.t_stop]
        tSpan = np.linspace(0, opts.t_stop)

        # #%% SOLUTION 1: Numerial solution of the model defined in reaction_kinetic_ODE_system using ode45
        #  sol1 = ode45(@(t,c) reaction_kinetic_ODE_system(c,k,params,n_particle,l), tSpan, ic, ode_opts)

        # reaction_kinetic_ODE_system(c,k,params,n_particle,l))
        # NOTE: dcdt_r is a wrapper for reaction_kinetic_ODE_system
        y0 = ic
        sol1 = scipy.integrate.solve_ivp(dcdt_r, (data_time[0], data_time[-1]), y0, method='RK45',
                                         # t_eval=tSpan,
                                         t_eval=data_time,
                                         # dense_output=False, events=None, vectorized=False,
                                         args=(k, params, n_particle, l),
                                         **ode_options)

        # #%% SOLUTION 2: Numerical solution of the diffusion model defined in diffusion_controlled_release
        # if l == 1:
        # % pure CuO dissolves without the diffusion controlled release
        #     sol2 = 0
        # else:
        # % Diffusion controlled release with final radius from solution
        # % 1 as starting point
        # r_final = sol1.y[4,:]
        r_final = sol1.y[4 - 1, -1]
        sol2 = diffusion_controlled_release(r_final, k, params, opts, c, l)

        # #%% Evaluation of solutions 1 and 2 at the experimental data (times)
        if 1:
            sim.append(sol1.y)
            # if l == 1-1:
            # % pure CuO dissolves without the diffusion controlled release
            #     sim2.append(0)  # Is this step still needed?
            #  else:
            sim2.append(sol2.y)
            # end

        # % Plot intermediate results for each iteration step
        # NOTE: plt.pause(0.01) can do this if needed, however, I would not recommend
        #      this in any simulation loop for performance considerations. - Stefan Endres
        if opts.disp_intermediate_results == 1:
            pass
            # Plots or save data:
            # figure(10)
            # semilogx(time_eval{l},sim{l}(2,:)+sim2{l}(end,:),'k-',time_eval{l},sim{l}(2,:),'k--',time_eval{l},sim2{l}(end,:),'k-.',time_eval{l},conc_eval{l},'s','MarkerEdgeColor',colors(l,:),'MarkerFaceColor',colors(l,:),'MarkerSize',10)
            # ylim([0 2.5])
            # h = legend('c_{Cu2+,sol1+sol2}','c_{Cu2+,sol1}','c_{Cu2+,sol2}','c_{Cu2+,exp}')
            # xlabel('Time t in h','FontSize',28)
            # ylabel('Concentration c_{Cu2+} in mM','FontSize',28)
            # set(findall(gcf,'Type','axes'),'LineWidth',2.835)
            # set(gca,'fontsize',30)
            # set(findall(gca, 'Type', 'Line'),'LineWidth',2.835)
            # set(gca,'XScale','log')
            # set(gca,'Xtick',[1 10 100])
            # xlim([0.5 opts.t_stop])
            # axis square
            # set(gcf,'Position',[50 50 650 650])
        else:
            pass

        # #%% Save model results for plot
        model['t'].append(sol1.t)
        # model.c_cu2p_1{l} = deval(sol1,sol1.x)
        model['c_cu2p_1'].append(sol1.y)
        # if l == 1-1:
        # model.c_cu2p_2{l} = 0
        #    model['c_cu2p_2'].append(0)
        # else:
        # model.c_cu2p_2{l} = deval(sol2,sol1.x)
        model['c_cu2p_2'].append(sol2.y)

        # #%% Calculation of the objective function, RMSE
        # minSquareError = minSquareError + sum((sim{l}(2,:) + sim2{l}(end,:) - conc_eval{l}).^2)/length(conc_eval{l})
        # NOTE: conc_eval became data_c_cu2p after the coversion
        minSquareError = np.linalg.norm(sol1.y - data_c_cu2p) / len(data_c_cu2p)

    return [minSquareError, model, sol2]
