%% Settings
function [opts, params, colors] = load_settings()
    % Options
    opts = struct("t_stop", 250,... % Max. time for evaluation in hrs
        "errorbars", 1,... % Display errorbars (0 - off, 1 - on)
        "N_nodes", 500,... % Number of nodes in numerical grid for solution of the diffusion equation (diffusion_controlled_release()) in -
        "t_max", 1e-6,... % Max. timestep for solution of diffusion equation in s
        "disp_intermediate_results", 1,... % Display intermediate results during itteration (0 - off, 1 - on)
        "disp_particle_profile", 1); % Display concentration profile in particles (c(r,t)) from diffusion equation (0 - off, 1 - on)
    % Note: N_nodes = 2500 required for sufficient accuracy
    
    % Physical constants and initial parameters
    params = struct("rho_CuO", 6.31*1e-21,... % Density of CuO in g/nm³
        "M_CuO", 79.545/1000,... % Molar mass of CuO in g/mmol
        "N_A", 6.022e23,... % Avogadro constant in #/mol
        "n_data", 4,... % Number of datasets for evaluation in -
        "c0", [2.5 0 5 0],... % Initial molar concentrations of [CuO, Cu2+, Amino Acid, Crystals] in mmol
        "d0", [11.8 12.0 10.3 10.7],... % Initial particle diameters as obtained from BET measurements in nm
        "m_r", 3.5,... % Reaction order in -
        "n_r", 4,... % Reaction order in -
        "ratioFe_Cu", [0 0.01 0.06 0.1]); % Iron-copper ratio in -
    
    % Color definitions
    black = [0 0 0]./255;
    cyan = [0 147 221]./255;
    green = [132 194 37]./255;
    red = [218 37 29]./255;
    colors(1,:) = black;
    colors(2,:) = red;
    colors(3,:) = cyan;
    colors(4,:) = green;
end