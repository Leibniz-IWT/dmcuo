%% Objective function that calculates the root mean square error (RMSE) between experimental data and model for a set of parameters (params_x)
%% The model data is obtained from the superimposed solution of the kinetically driven relase (solution 1) and the diffusion limited release (solution 2)
function [minSquareError, model, sol2] = obj_fun(params_x,data,params,opts,colors)

    %% Error initialization
    minSquareError = 0.0;
    
    %% Model calculation for pure CuO, 1%Fe+CuO, 6%Fe+CuO and 10%Fe+CuO
    for l = 1:params.n_data

        %% Prepare experimental data - Concentration profiles
        % Indicies of the data (table with varying length) 
        first_idx = 1;
        last_idx = find(data.time(l,:) ~= 0,1,'last');
        % Remove empty datapoints (zeros) in the table
        time_eval{l} = data.time(l,first_idx:last_idx);
        conc_eval{l} = data.c_cu2p(l,first_idx:last_idx);
        % Remove datapoints exceeding opts.t_stop limit
        time_eval{l} = time_eval{l}(find(time_eval{l} <= opts.t_stop)); 
        conc_eval{l} = conc_eval{l}(find(time_eval{l} <= opts.t_stop));
        
        %% Initial conditions and parameters for the ODE system describing the kinetics (ode45 solver)
        % Initial particle diameter as obtained from BET measurements
        d0 = params.d0(l); % in nm
        % Initial particle radius
        r0 = 0.5*d0; % in nm
        
        %% Fit parameters for the ODE system
        % Dissolution rate constants
        k(1) = params_x(4)*1e-36;
        k(2) = params_x(5)*1e-36;
        k(3) = params_x(6)*1e-36;
        k(4) = params_x(7)*1e-36;
        % Diffusion coefficient
        k(5) = params_x(2);
        % Atomic surface concentration per unit surface area in #/nm?
        k(6) = params_x(1)*1000*params.M_CuO/params.rho_CuO/params.N_A/0.073;
        
        %% Initial concentrations 
        c(1) = params_x(3); % C_CuO (variable) in mM
        c(2) = params.c0(2); % C_Cu2+ (fixed) in mM
        c(3) = params.c0(3); % C_AA (amino acid, fixed) in mM
        
        %% Calculated values based on input
        V0 = c(1)*params.M_CuO/params.rho_CuO; % Particle volume in nm?
        n_particle = V0/(4/3*pi()*r0*r0*r0); % Number of particles in -
        a0 = 3*V0/r0; % Total particle surface area in nm?
        c(4) = r0; % r_CuO - Particle radius in nm
        c(5) = a0; % a_CuO - Particle surface area in nm?
        c(6) = V0; % V_CuO - Particle volume in nm?
        
        %% Initial iron-copper and copper-iron ratio
        c(7) = params.ratioFe_Cu(l)*c(5)/k(6); % f_fe
        c(8) = (1-params.ratioFe_Cu(l))*c(5)/k(6); % f_Cu
        
        %% Set initial conditions for the ODE system
        ic = [c(1) c(2) c(3) c(4) c(5) c(6) c(7) c(8)];

        %% Options for ode45
        ode_opts = odeset('RelTol',2e-14,'AbsTol',1e-12,'Nonnegative',[]);

        %% Timespan of model simulation
        tSpan = [0:1:opts.t_stop];
        
        %% SOLUTION 1: Numerial solution of the model defined in reaction_kinetic_ODE_system using ode45
        sol1 = ode45(@(t,c) reaction_kinetic_ODE_system(c,k,params,n_particle,l), tSpan, ic, ode_opts);
        
        %% SOLUTION 2: Numerical solution of the diffusion model defined in diffusion_controlled_release
        if l == 1 
            % pure CuO dissolves without the diffusion controlled release
            sol2 = 0;
        else
            % Diffusion controlled release with final radius from solution
            % 1 as starting point
            r_final = sol1.y(4,end);
            sol2 = diffusion_controlled_release(r_final,k,params,opts,c,l);
        end
        
        %% Evaluation of solutions 1 and 2 at the experimental data (times)
        sim{l} = deval(sol1,time_eval{l});
        if l == 1
            % pure CuO dissolves without the diffusion controlled release
            sim2{l} = 0;
        else
            sim2{l} = deval(sol2,time_eval{l});
        end
        
        % Plot intermediate results for each itteration step
        if opts.disp_intermediate_results == 1
            figure(10)
            semilogx(time_eval{l},sim{l}(2,:)+sim2{l}(end,:),'k-',time_eval{l},sim{l}(2,:),'k--',time_eval{l},sim2{l}(end,:),'k-.',time_eval{l},conc_eval{l},'s','MarkerEdgeColor',colors(l,:),'MarkerFaceColor',colors(l,:),'MarkerSize',10)
            ylim([0 2.5])
            h = legend('c_{Cu2+,sol1+sol2}','c_{Cu2+,sol1}','c_{Cu2+,sol2}','c_{Cu2+,exp}');
            xlabel('Time t in h','FontSize',28);
            ylabel('Concentration c_{Cu2+} in mM','FontSize',28);
            set(findall(gcf,'Type','axes'),'LineWidth',2.835);
            set(gca,'fontsize',30);
            set(findall(gca, 'Type', 'Line'),'LineWidth',2.835);
            set(gca,'XScale','log');
            set(gca,'Xtick',[1 10 100]);
            xlim([0.5 opts.t_stop]);
            axis square;
            set(gcf,'Position',[50 50 650 650])
        else
        end
        
        %% Save model results for plot
        model.t{l} = sol1.x;
        model.c_cu2p_1{l} = deval(sol1,sol1.x);
        if l == 1
            model.c_cu2p_2{l} = 0;
        else
            model.c_cu2p_2{l} = deval(sol2,sol1.x);
        end
    
        %% Calculation of the objective function, RMSE
        minSquareError = minSquareError + sum((sim{l}(2,:) + sim2{l}(end,:) - conc_eval{l}).^2)/length(conc_eval{l});
    end
end
