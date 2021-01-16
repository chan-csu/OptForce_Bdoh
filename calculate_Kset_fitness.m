function y = calculate_Kset_fitness(Experimental_Data, solutions, ...
    screening_conditions, S)

    %%% function calculates fitness score for a parameter set flux
    %%% prediction compared to experimentally observed fluxes (stored in
    %%% Experimental_Data
    
    ref_fluxes = Experimental_Data.ref_fluxes;
    flux_rxns = Experimental_Data.flux_rxns;
    conditions = Experimental_Data.conditions;
    weighted_rxn = Experimental_Data.weighted_rxn;
    
    y_iter = 0;
                
    for m = 1:length(screening_conditions)
        perturbed_fluxes = solutions(:,m);
        
        %%% penalize fitness score if parameter set prediction is not a
        %%% valid mass-balanced steady state solution (ex: integration did
        %%% not complete over specified timeframe and new steady state 
        %%% not-reached)
        if max(abs(S*perturbed_fluxes))>= 0.1
            y_iter = 100 + y_iter;     
        
        %%%% penalize fitness score if no flux solution predicted for
        %%%% parameter set (i.e. failed integration)
        elseif sum(perturbed_fluxes) == 0
            y_iter = 100 + y_iter;
                
        else
            
            predicted_flux = perturbed_fluxes(flux_rxns(3:end));                   
                   
            ref_flux = ref_fluxes(3:end,conditions == ...
                screening_conditions(m));
            weighted_rxn_flux = abs(ref_flux(flux_rxns == weighted_rxn));
            n_refs = length(ref_flux);
        
            y_iter_next = 1/n_refs*...
                sum((abs(predicted_flux - ref_flux)./weighted_rxn_flux)); 
            y_iter = y_iter + y_iter_next;
                
        end
           
        y =  1/length(screening_conditions) * y_iter;
        
    end   
    
end