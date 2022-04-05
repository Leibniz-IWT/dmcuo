import numpy as np
import scipy.integrate
from radial_diffusion_PDE import radial_diffusion_PDE

#% Function to calculate diffusion limited release
def diffusion_controlled_release(r_f, k, params, opts, c, l):
    #% Initial CuO concentration
    #c0 = c(1)
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
  #  y_0(1:N+1, 1) = c_S0 #% in #/nm�
  #  y_0(N+1, 1) = 0 #% in #/nm�
  #  y_0(N+2, 1) = r_f #% R(t) in nm
  #  y_0(N+3, 1) = 0 #% c_cu2+ in mM
    y_0 = np.ones(N + 3)
    y_0[:N] = c_S0
    y_0[N] = 0
    y_0[N + 1] = r_f
    y_0[N + 2] = 0

    #% Time for evaluation in h
  #  tSpan = [0:1:opts.t_stop]
    #tSpan = np.linspace(0, opts.t_stop, int(opts.t_stop) + 1)
    t_span = (0, opts.t_stop)
    # from : x = [a:b: c] and y = linspace(a, c, (c - a) / b + 1)?
    
    #% ODE solver, stiff ode15s solver chosen due to strong concentration
    #% gradients in the beginning
  #  ode_opts = odeset('RelTol',2e - 10,'AbsTol',1e - 12,'Nonnegative',[],'InitialStep',1e - 20)
  #  sol2 = ode15s(@(t,y) radial_diffusion_PDE(t,y,params,opts,k,l,c0), tSpan, y_0, ode_opts)

    args = (params, opts, k, l, c0)
    # first_step=None, min_step=0.0, max_step=inf, rtol=0.001, atol=1e-06, jac=None, lband=None, uband=None,
    options = {'first_step': 1e-20,
              'atol': 1e-12,
               # 'atol': 1e-3,
               'rtol': 1e-10,
               #'rtol': 1e-2
               }
    #NOTE: tSpan can be added to t_eval here for interpolations at those points:
    sol2 =scipy.integrate.solve_ivp(radial_diffusion_PDE, t_span, y_0, method='LSODA',
                                   # t_eval=None, dense_output=False, events=None, vectorized=False,
                                    args=args, **options)

    # fun(t, y). Here t is a scalar, and there are two options for the ndarray y
    return sol2