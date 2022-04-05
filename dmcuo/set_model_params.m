function [params_0, params_lb, params_ub] = set_model_params(k1,k2,k3,k4,ks,D_solid,c1,c2,c3,c4)
    % Parameter definition and upper and lower boundaries
    params_0(1) = ks; params_lb(1) = 1.0; params_ub(1) = 1.75;
    params_0(2) = D_solid; params_lb(2) = D_solid/5; params_ub(2) = D_solid*5;
    params_0(3) = c1; params_lb(3) = 2.0; params_ub(3) = 2.25;
    params_0(4) = k1; params_lb(4) = k1/2; params_ub(4) = k1*4;
    params_0(5) = k2; params_lb(5) = k2/2; params_ub(5) = k2*4;
    params_0(6) = k3; params_lb(6) = k3/2; params_ub(6) = k3*4;
    params_0(7) = k4; params_lb(7) = k4/2; params_ub(7) = k4*4;
    
    % Parameter fitting runs for approximately 2 days if N = 2500 nodes is chosen. Hence, final values
    % obtained from the fitting procedure are included
    params_0(1) = 1.106424465; params_lb(1) = 1.106424465; params_ub(1) = 1.106424465;
    params_0(2) = 0.000019913369467558; params_lb(2) = 0.000019913369467558; params_ub(2) = 0.000019913369467558;
    params_0(3) = 2.04372584529584; params_lb(3) = 2.04372584529584; params_ub(3) = 2.04372584529584;
    params_0(4) = 20.0001185922097; params_lb(4) = 20.0001185922097; params_ub(4) = 20.0001185922097;
    params_0(5) = 4.14276079865522; params_lb(5) = 4.14276079865522; params_ub(5) = 4.14276079865522;
    params_0(6) = 2.85563976205523; params_lb(6) = 2.85563976205523; params_ub(6) = 2.85563976205523;
    params_0(7) = 0.678991631490014; params_lb(7) = 0.678991631490014; params_ub(7) = 0.678991631490014; 
end
