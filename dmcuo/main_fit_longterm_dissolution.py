#!/usr/bin/env python3
"""
## Main routine for fitting release kinetic model to experimental concentration profiles in Figure 2.20 (A) and (B)
## Author: Hendrik Naatz
"""
# Library imports:
import numpy as np
import codecs, json
# Non std.:
import scipy.io
import scipy.optimize
import matplotlib.pyplot as plt

# Module specific imports:
from load_settings import load_settings
from set_model_params import set_model_params
from set_model_params import set_model_params
#from obj_fun import obj_fun
from reaction_kinetic_ODE_system import reaction_kinetic_ODE_system
from diffusion_controlled_release import diffusion_controlled_release
data_path = './Data/dissolution_profiles_cu2p.json'

## Procedure
# 1. Plot experimental data, i.e. concentration files
# 2. Fit of kinetic model
# 3. Plot experimental data together with model fit

## Clean workspace and cmd window
#clear all clc warning off

## Load settings
# Note: N_Nodes should be 2500 for sufficient precision
opts, params, colors = load_settings()

## Read in dissolution profiles, i.e. sampling times (time) in hours
# and dissolved Cu2+ species (c_cu2p) in mM, obtained via UV-Vis
# measurements according to [DOI: 10.1002/anie.201916183]
# load data from json:
obj_text = codecs.open(data_path, 'r', encoding='utf-8').read()
data = json.loads(obj_text)
print(f'data.keys() = {data.keys()}')
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

#  clean the data of the concentration profiles [new for Python translation -SCE]:
#data['mean']
#data.time

# Test function conversions:
if 0:
    # Test ODE system:
    l =0
    k = np.ones(6)
    n_particle = 1e10
    c = np.ones(8)
    print(reaction_kinetic_ODE_system(c,k,params,n_particle,l))
    # Test diffusion function
    r_f = 1  # r_final, final concentration?
    print(diffusion_controlled_release(r_f, k, params, opts, c, l))
    # obj_fun(params_x,data,params,opts,colors)
# % (In)equalitiy condition
Aeq = []
Beq = []

# % Definition of fit function (RMSE) to minimize in file obj_fun
#  objective = @(params_x)obj_fun(params_x,data,params,opts,colors)
# % Fitting of model parameters with minimum root mean square error (params_f) using
# % fmincon
#[params_f, Fend] = fmincon(objective, params_0, Aeq, Beq, [], [], params_lb, params_ub)

#scipy.optimize.least_squares

# % Calculation of model profiles using final params_f obtained from fmincon
#[minSquareError, model, sol2] = obj_fun(params_f, data, params, opts, colors)

### Plots
## Plot experimental data (concentration profiles)
fig = plt.figure(1)
print(data['time'][0].shape)
markers = ['o', '^', 's', 'h']
ax = plt.gca()
for i in range(params.n_data):
    print(i)
    t = data['time'][i]
    y = data['mean'][i]
    ax.scatter(t, y, marker=markers[i], s=30, color=colors[i],
                label=data_keys[i])
    yerror = data['max'][i] - data['min'][i]
    ax.errorbar(t, y, yerr=yerror
                 , ls='none', color='Black',
                 elinewidth=2, capthick=2, errorevery=1, alpha=1, ms=4, capsize=5
                 )
ax.set_xscale('log')
ax.set_xlabel('Time t in h'#, fontsize=28
              )
ax.set_ylabel(r'Concentration $c_{Cu2+}$ in mM'#, fontsize=28
              )
plt.legend()

# %#% Plot final model fit and experimental data (concentration profiles)
if 0:
    figure(2)
    for i in range(1,params.n_data):
       # plot(model.t{i},model.c_cu2p_1{i}(2,:)+model.c_cu2p_2{i}(end,:),'k-')
       # hold on
       # e = errorbar(data.time(i,:),data.c_cu2p(i,:),opts.errorbars*(data.c_cu2p_max(i,:)-data.c_cu2p_min(i,:)),'s','MarkerEdgeColor',colors(i,:),'MarkerFaceColor',colors(i,:),'MarkerSize',10)
        e.CapSize = 12
        e.LineWidth = 2
        e.Color = 'black'

    h = legend('','CuO','','CuO+1#%Fe','','CuO+6#%Fe','','CuO+10#%Fe')
    xlabel('Time t in h','FontSize',28)
    ylabel('Concentration c_{Cu2+} in mM','FontSize',28)
    set(findall(gcf,'Type','axes'),'LineWidth',2.835)
    set(gca,'fontsize',30)
    set(findall(gca, 'Type', 'Line'),'LineWidth',2.835)
    set(gca,'XScale','log')
    set(gca,'Xtick',[1, 10, 100])
    xlim([0.5, opts.t_stop])
   # axis square
    set(gcf,'Position',[50, 50, 650, 650])


    #%#% Plot copper concentration profile in particle c_cu(r(t),t)
    if opts.disp_particle_profile == 1:
        for i in range( 1,length(sol2.x)):
            figure(3)
            radius = sol2.y(end-1,i) #% Shrinking particle size, i.e. moving boundary condition, R(t)
          #  c_cu_r = sol2.y(1:end-2,i) #% copper concentration profile in particle c_cu(r,t)
          #  plot([1:1:length(sol2.y)-2]/(length(sol2.y)-2)*radius,sol2.y(1:end-2,i)/max(sol2.y(1:end-2,i))) #% c_cu(R(t)), normalized
          #  text(0.5,0.1, ['t = ' num2str(sol2.x(i)) ' h'],'FontSize',28)
            xlabel('Radius r in nm','FontSize',28)
            ylabel('Norm. concentration c_{Cu}(R(t),t)','FontSize',28)
            set(findall(gcf,'Type','axes'),'LineWidth',2.835)
            set(gca,'fontsize',30)
            set(findall(gca, 'Type', 'Line'),'LineWidth',2.835)
            xlim([0, 5])
           # axis square
            pause(0.01)
            set(gcf,'Position',[50, 50, 650, 650])


    disp("Fit Parameter:")
    params_f

plt.show()