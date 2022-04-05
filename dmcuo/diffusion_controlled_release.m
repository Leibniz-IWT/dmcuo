% Function to calculate diffusion limited release
function sol2 = diffusion_controlled_release(r_f,k,params,opts,c,l)
    % Initial CuO concentration 
    c0 = c(1);
    
    % Spacing of the numberical grid
    N = opts.N_nodes; % in - 
    
    % Initial profile y_0 for method of lines
        % 1:N - concentration,
        % N+1 - surface concentration
        % N+2 - radius
        % N+3 - c_cu2+ concentration
    f_cu_fe = (1-params.ratioFe_Cu(l)); % Cu/(Cu+Fe) atomic ratio in -
    c_p = params.rho_CuO/params.M_CuO*params.N_A; % in #/nm³
    c_S0 = c_p*f_cu_fe; % in #/nm³
    y_0(1:N+1,1) = c_S0; % in #/nm³
    y_0(N+1,1) = 0; % in #/nm³
    y_0(N+2,1) = r_f; % R(t) in nm
    y_0(N+3,1) = 0; % c_cu2+ in mM
    
    % Time for evaluation in h
    tSpan = [0:1:opts.t_stop];
    
    % ODE solver, stiff ode15s solver chosen due to strong concentration
    % gradients in the beginning
    ode_opts = odeset('RelTol',2e-10,'AbsTol',1e-12,'Nonnegative',[],'InitialStep',1e-20);
    sol2 = ode15s(@(t,y) radial_diffusion_PDE(t,y,params,opts,k,l,c0), tSpan, y_0, ode_opts);
end