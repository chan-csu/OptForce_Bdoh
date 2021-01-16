function norm_slope = calc_norm_slope(X,K,S_f_b,uptake_rxns,uptake_values)
    %%% function calculates norm of all dx/dt values for a given parameter
    %%% set at a specific metabolite concentration profile; this function
    %%% allows us to determine if a steady state has been reached during
    %%% integration as we can see when all dx/dt values approach 0
    
    %%% calculate elementary reaction fluxes

    v = K;

    elem_rxns = 1:length(v);
    elem_rxns(uptake_rxns) = [];

    for j = 1:length(elem_rxns)

        for i=1:length(X) 
            if S_f_b(i,elem_rxns(j)) < 0
                v(elem_rxns(j))=v(elem_rxns(j))*...
                    (X(i)^abs(S_f_b(i,elem_rxns(j)))); 
            end
        end
    
    end
    
    %%% DAE implementation: set uptake reaction fluxes to experimentally
    %%% observed values
    v(uptake_rxns) = uptake_values;


%%%%%%%%%%%%%%%%%%
dx = S_f_b*v;

norm_slope = norm(dx);