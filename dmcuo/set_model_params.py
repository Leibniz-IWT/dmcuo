import numpy as np

def set_model_params(k1, k2, k3, k4, ks, D_solid, c1, c2, c3, c4):
    # Parameter definition and upper and lower boundaries
    # (NOTE: In converting the .m code it was faster to macro the indices tp minus one
    params_0 = np.zeros(7)
    params_lb = np.zeros(7)
    params_ub = np.zeros(7)
    params_0[1-1] = ks
    params_lb[1-1] = 1.0
    params_ub[1-1] = 1.75
    params_0[2-1] = D_solid
    params_lb[2-1] = D_solid / 5
    params_ub[2-1] = D_solid * 5
    params_0[3-1] = c1
    params_lb[3-1] = 2.0
    params_ub[3-1] = 2.25
    params_0[4-1] = k1
    params_lb[4-1] = k1 / 2
    params_ub[4-1] = k1 * 4
    params_0[5-1] = k2
    params_lb[5-1] = k2 / 2
    params_ub[5-1] = k2 * 4
    params_0[6-1] = k3
    params_lb[6-1] = k3 / 2
    params_ub[6-1] = k3 * 4
    params_0[7-1] = k4
    params_lb[7-1] = k4 / 2
    params_ub[7-1] = k4 * 4

    # Parameter fitting runs for approximately 2 days if N = 2500 nodes is chosen. Hence, final values
    # obtained from the fitting procedure are included
    params_0[1-1] = 1.106424465
    params_lb[1-1] = 1.106424465
    params_ub[1-1] = 1.106424465
    params_0[2-1] = 0.000019913369467558
    params_lb[2-1] = 0.000019913369467558
    params_ub[2-1] = 0.000019913369467558
    params_0[3-1] = 2.04372584529584
    params_lb[3-1] = 2.04372584529584
    params_ub[3-1] = 2.04372584529584
    params_0[4-1] = 20.0001185922097
    params_lb[4-1] = 20.0001185922097
    params_ub[4-1] = 20.0001185922097
    params_0[5-1] = 4.14276079865522
    params_lb[5-1] = 4.14276079865522
    params_ub[5-1] = 4.14276079865522
    params_0[6-1] = 2.85563976205523
    params_lb[6-1] = 2.85563976205523
    params_ub[6-1] = 2.85563976205523
    params_0[7-1] = 0.678991631490014
    params_lb[7-1] = 0.678991631490014
    params_ub[7-1] = 0.678991631490014
    return params_0, params_lb, params_ub
