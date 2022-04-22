#!/usr/bin/env python3
"""
## Main routine for fitting release kinetic model to experimental concentration profiles in Figure 2.20 (A) and (B)
## Author: Hendrik Naatz
"""
# Library imports:
import numpy as np
import codecs, json
import copy
# Non std.:
import scipy.io
import scipy.optimize
import matplotlib.pyplot as plt

# Module specific imports:
from dmcuo.load_settings import load_settings
from dmcuo.set_model_params import set_model_params
#from obj_fun import obj_fun
from dmcuo.reaction_kinetic_ODE_system import reaction_kinetic_ODE_system
from dmcuo.diffusion_controlled_release import diffusion_controlled_release
from dmcuo.obj_fun import obj_fun
data_path = 'dmcuo/Data/dissolution_profiles_cu2p.json'

## Procedure
# 1. Plot experimental data, i.e. concentration files
# 2. Fit of kinetic model
# 3. Plot experimental data together with model fit

## Load settings
# Note: N_Nodes should be 2500 for sufficient precision
opts, params, colors = load_settings()

## Read in dissolution profiles, i.e. sampling times (time) in hours
# and dissolved Cu2+ species (c_cu2p) in mM, obtained via UV-Vis
# measurements according to [DOI: 10.1002/anie.201916183]
# load data from json:
obj_text = codecs.open(data_path, 'r', encoding='utf-8').read()
data = json.loads(obj_text)
for d in data:
    data[d] = np.array(data[d])  # convert lists back to numpy arrays

#  Keys if the 4 data sets for labels
data_keys = ['CuO', 'CuO + 1% Fe', 'CuO + 6% Fe', 'CuO + 10% Fe']
## Fit of kinetic model
# Initial guess of model parameters
k1 = 20  #% Dissolution rate constant for 0#% Fe
k2 = 4  #% Dissolution rate constant for 1#% Fe
k3 = 2.75  #% Dissolution rate constant for 6#% Fe
k4 = 1  #% Dissolution rate constant for 10#% Fe
ks = 1.25  #% Atomic surface concentration
D_solid = 0.00005  #% Diffusion coefficient
c1 = 2.25  # % C_CuO,0 inital concentration, always close to inital concentration given

#% Initial values
c2 = 0  #% C_Cu2+,0 initial concentration
c3 = 5  #% C_AminoAcid,0 initial concentration
c4 = 0  #% C_crystalization,0 initial concentration

# % Set initial model parameters and boundaries (currently fixed to optimized
# % values)
params_0, params_lb, params_ub = set_model_params(k1, k2, k3, k4, ks, D_solid, c1, c2, c3, c4)
# % Number of fit parameters as defined in params_0
nvar = len(params_0)
# % Definition of fit function (RMSE) to minimize in file obj_fun
#  objective = @(params_x)obj_fun(params_x,data,params,opts,colors)
obj_fun(params_0, data, params, opts, colors)
# % Fitting of model parameters with minimum root mean square error (params_f)
bounds = []
for plb, pub in zip(params_lb, params_ub):
    bounds.append([plb, pub])
if 0:
    def obj(x,data, params, opts, colors):
        params_0 = x
        obj = obj_fun(params_0, data, params, opts, colors)
        return obj[0]

    res = scipy.optimize.minimize(obj, params_0, args=(data, params, opts, colors))
    #scipy.optimize.least_squares
    params_f = res.x
else:
    params_f = params_0

# % Calculation of model profiles using final params_f obtained from fmincon
tSpan = np.linspace(0, data['time'][-1][-1], 3000)
tl = data['time'].shape[0]
data_time = copy.copy(data['time'])
data['time'] = []
for l in range(tl):
    data['time'].append(tSpan)

[minSquareError, model, sol2] = obj_fun(params_f, data, params, opts, colors)

data['time'] = data_time
### Plots
## Plot experimental data (concentration profiles)
fig = plt.figure(1)
markers = ['o', '^', 's', 'h']
ax = plt.gca()
for i in range(params.n_data):
    t = data['time'][i]
    y = data['mean'][i]
    ax.scatter(t, y, marker=markers[i], s=30, color=colors[i],
                label=data_keys[i])
    yerror = data['max'][i] - data['min'][i]
    ax.errorbar(t, y, yerr=yerror, ls='none', color='Black',
                 elinewidth=2, capthick=2, errorevery=1, alpha=1, ms=4, capsize=5
                 )
ax.set_xscale('log')
ax.set_xlabel('Time t in h')
ax.set_ylabel(r'Concentration $c_{Cu2+}$ in mM')
plt.legend()

# %#% Plot final model fit and experimental data (concentration profiles)
fig = plt.figure(2)
markers = ['o', '^', 's', 'h']
ax = plt.gca()
for i in range(params.n_data):
    print(i)
    t = data['time'][i]
    y = data['mean'][i]
    ax.scatter(t, y, marker=markers[i], s=30, color=colors[i],
               label=data_keys[i])
    yerror = data['max'][i] - data['min'][i]
    ax.errorbar(t, y, yerr=yerror, ls='none', color='Black',
                elinewidth=2, capthick=2, errorevery=1, alpha=1, ms=4, capsize=5
                )
    ax.plot(model['t'][i][1:], model['c_cu2p_1'][i][2-1][1:], '--', color=colors[i])
    ax.plot(model['t'][i][1:], model['c_cu2p_1'][i][2-1][1:] + model['c_cu2p_2'][i][1:], color=colors[i])
ax.set_xscale('log')
ax.set_xlabel('Time t in h'  # , fontsize=28
              )
ax.set_ylabel(r'Concentration $c_{Cu2+}$ in mM'  # , fontsize=28
              )
plt.legend()

print(f'Final parameter set:')
print(f'params_f = {params_f}')

plt.show()