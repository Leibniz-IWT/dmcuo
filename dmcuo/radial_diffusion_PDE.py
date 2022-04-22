import numpy as np

#%#% Calcuation of the diffusion limited release
def radial_diffusion_PDE(t ,y , params, opts, k, l, c0):
    #% Parameter definition
    # Translation note: all c(_) indices where changed to c[_ - 1] t keep track of originals with macros - SCE
    r0 = params.d0[l]/2  #% in nm
    D_solid = k[5 - 1]  #% in nm²/h

    # Add in r0 calculation:
    d0 = params.d0[l]  # % in nm
    r0 = 0.5 * d0  # % Initial particle radius in nm

    #% Parameter
    N = len(y) - 3
    dz = 1/N
    c = y  # c[1:N+3) = y;

    #% Inititialize dcdt
    dcdt = np.zeros(N + 3)

    #% Symmetry BC
    R_t = c[N + 2 - 1]  # R(t) = y[N + 1] = c[N + 1]
    dcdt[0] = 6 * D_solid/(R_t**2) * ((c[2 - 1] - c[1 - 1])/(dz**2))

    #% Spherical particle, Discretization of PDE using method of lines (MOL)
    for i in range(2, N):  # All elements between the second (1) and N interior element
        dRdt = D_solid * (params.M_CuO/(params.rho_CuO * params.N_A)) * (
                (3 * c[N + 1 - 1] - 4 * c[N - 1] + c[N - 1 - 1])/(2 * R_t * dz))
        dcdt[i] = ((((i - 1)/(2 * R_t)) * dRdt + D_solid / ((i - 1) * (dz * R_t)**2))
                    * (c[i + 1 - 1] - c[i - 1 - 1])
                   + D_solid / ((dz * R_t) **2) * (c[i + 1 - 1] - 2 * c[i - 1] + c[i - 1 - 1])
                  )  # % Check
    dcdt = np.nan_to_num(dcdt)
    #% Surface condition
    dcdt[N] = 0

    #% Moving interface BC
    dcdt[N + 1] = D_solid * params.M_CuO/params.rho_CuO/params.N_A/c[N + 2 - 1] * (3 * 0
                  - 4 * c[N- 1] + c[N - 1 - 1])/2/dz

    R_t = c[N + 2 - 1]
    dcdt[N + 2] = - 4 * np.pi * R_t * R_t * D_solid * params.M_CuO/params.rho_CuO/(params.N_A * R_t) * (3 * 0
                  - 4 * c[N - 1] + c[N - 1 - 1])/(2*dz) * params.rho_CuO/params.M_CuO * c0 * params.M_CuO/params.rho_CuO/(4/3
                  * np.pi * r0 * r0 * r0) * 0.5

    return dcdt