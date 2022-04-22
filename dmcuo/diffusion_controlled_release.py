import numpy as np
import scipy.integrate
from dmcuo.radial_diffusion_PDE import radial_diffusion_PDE

#% Function to calculate diffusion limited release
def diffusion_controlled_release(r_f, k, params, opts, c, l, t_eval=None):
    #% Initial CuO concentration
    c0 = c[0]

    #% Spacing of the numberical grid
    N = opts.N_nodes #% in  -
    
    #% Initial profile y_0 for method of lines
        #% 1:N  -  concentration,
        #% N+1  -  surface concentration
        #% N+2  -  radius
        #% N+3  -  c_cu2+ concentration
    f_cu_fe = (1 - params.ratioFe_Cu[l]) #% Cu/(Cu+Fe) atomic ratio in  -
    c_p = params.rho_CuO/params.M_CuO * params.N_A #% in #/nm�
    c_S0 = c_p * f_cu_fe #% in #/nm�
    y_0 = np.ones(N + 3)
    y_0[:N] = c_S0  # in #/nm�
    y_0[N] = 0  # in #/nm�
    y_0[N + 1] = r_f  # R(t) in nm
    y_0[N + 2] = 0  #  c_cu2+ in mM

    #% Time for evaluation in h
    t_span = (0, opts.t_stop)
    
    #% ODE solver, stiff ode15s solver chosen due to strong concentration
    #% gradients in the beginning


    args = (params, opts, k, l, c0)
    # first_step=None, min_step=0.0, max_step=inf, rtol=0.001, atol=1e-06, jac=None, lband=None, uband=None,
    options = {'first_step': 1e-20,
               'atol': 1e-12,
               #'atol': 1e-3,
               'rtol': 1e-10,
              # 'rtol': 1e-2
               }

    sol2 = scipy.integrate.solve_ivp(radial_diffusion_PDE, t_span, y_0,
                                     #method='LSODA',
                                       method='RK45',
                                    t_eval=t_eval, dense_output=False, events=None, vectorized=False,
                                    args=args,
                                    # **options
                                     )

    # fun(t, y). Here t is a scalar, and there are two options for the ndarray y
    return sol2