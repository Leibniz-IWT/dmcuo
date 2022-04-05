%% Calcuation of the diffusion limited release
function fval = radial_diffusion_PDE(t,y,params,opts,k,l,c0)

% Parameter definition
r0 = params.d0(l)/2; % in nm
D_solid = k(5); % in nm²/h

% Parameter
N = length(y)-3;
dz = 1/(N); 
c(1:N+3) = y;

% Inititialize dcdt
dcdt = zeros(N+3,1);

% Symmetry BC
dcdt(1) = 6*D_solid/c(N+2)/c(N+2)/dz/dz*(c(2)-c(1));

% Spherical particle, Discretization of PDE using method of lines (MOL)
for i=2:N
    dcdt(i) = ((i-1)/2/c(N+2)*(D_solid*params.M_CuO/params.rho_CuO/params.N_A/c(N+2)*(3*0-4*c(N)+c(N-1))/2/dz)+D_solid/(i-1)/dz/dz/c(N+2)/c(N+2))*(c(i+1)-c(i-1))+D_solid/dz/dz/c(N+2)/c(N+2)*(c(i+1)-2*c(i)+c(i-1)); % Check
end

% Surface condition
dcdt(N+1) = 0;

% Moving interface BC
dcdt(N+2) = D_solid*params.M_CuO/params.rho_CuO/params.N_A/c(N+2)*(3*0-4*c(N)+c(N-1))/2/dz;

% dc_cu2+dt
dcdt(N+3) = -4*pi()*c(N+2)*c(N+2)*D_solid*params.M_CuO/params.rho_CuO/params.N_A/c(N+2)*(3*0-4*c(N)+c(N-1))/2/dz*params.rho_CuO/params.M_CuO*c0*params.M_CuO/params.rho_CuO/(4/3*pi()*r0*r0*r0);

% Extract fval from dc/dt
fval = dcdt(1:N+3);
end