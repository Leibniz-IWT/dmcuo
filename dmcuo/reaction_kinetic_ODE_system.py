import numpy as np
#%#% Reaction kinetic ODE system calculates concnetration profiles and particle diameter, surface area and size
def reaction_kinetic_ODE_system(c, k, params, n_particle, l):
    #% Partial reaction orders
    m = params.m_r
    n = params.n_r
    
    #% Reaction kinetic ODE system
    dcdt = np.zeros(8)
    #%(1) C_CuO
    dcdt[1-1] = -k[l] * (1/k[5] * (c[8-1]/(c[7-1] + c[8-1])) * c[5-1]) ** (m/2) * c[3-1] ** (n/2)
    #%(2) C_Cu2 +
    dcdt[2-1] = k[l] * (1/k[5] * (c[8-1]/(c[7-1] + c[8-1])) * c[5-1]) ** (m/2) * c[3-1] ** (n/2)
    #%(3) C_Amino acid
    dcdt[3-1] = -2 * k[l] * (1/k[5] * (c[8-1]/(c[7-1] + c[8-1])) * c[5-1]) ** (m/2) * c[3-1] ** (n/2)
    #%(4) r_CuO
    dcdt[4-1] = -k[l]/c[5-1] * (1/k[5] * (c[8-1]/(c[7-1] + c[8-1])) * c[5-1]) ** (m/2)  \
                   * c[3-1] ** (n/2) * params.M_CuO/params.rho_CuO
    #%(5) a_CuO
    dcdt[5-1] = -n_particle/c[5-1] * 8 * np.pi * c[4-1] * k[l] * (1/k[5] * (c[8-1]/(c[7-1] + c[8-1]))  \
                   * c[5-1]) ** (m/2) * c[3-1] ** (n/2) * params.M_CuO/params.rho_CuO
    #%(6) V_CuO
    dcdt[6-1] = 0 #% Volume isn't required for rest of calculation
    #%(6) n_Fe
    dcdt[7-1] = k[l] * (1/k[5] * (c[8-1]/(c[7-1] + c[8-1])) * c[5-1]) ** (m/2) * c[3-1] ** (n/2)  \
                   * params.ratioFe_Cu[l]/(1-params.ratioFe_Cu[l]) * params.N_A/1000
    #%(7) n_cu
    dcdt[8-1] = -1/k[5]/c[5-1] * n_particle * 8 * np.pi * c[4-1] * k[l] * (1/k[5] * (c[8-1]/(c[7-1] + c[8-1])) \
                   * c[5-1]) ** (m/2) * c[3-1] ** (n/2) * params.M_CuO/params.rho_CuO-k[l] * (1/k[5] \
                   * (c[8-1]/(c[7-1] + c[8-1])) * c[5-1]) ** (m/2) * c[3-1] ** (n/2) * params.ratioFe_Cu[l] \
                   /(1-params.ratioFe_Cu[l]) * params.N_A/1000
    return dcdt