function [fvals, all_solutions, all_complete_times, ...
    all_ode_warn_flags, all_slope_norms] = ... 
    calculate_perturbed_steady_states(Network_Data, Experimental_Data, ...
    screening_conditions, All_K, All_fractions, t_interval, uptake_rxns)

    fprintf(1,'\n');
    fprintf(1,'Calculate perturbed steady state flux results for ensemble of k-sets...');
    fprintf(1,'\n');  
    
    nPerturbs = length(screening_conditions);
    nKsets = size(All_K,2);
        
    fvals = zeros(nKsets,1);
    solutions = zeros(size(Network_Data.S,2),nPerturbs);
    all_solutions = cell(nPerturbs,1);
    all_complete_times = zeros(nKsets, nPerturbs);
    all_ode_warn_flags = zeros(nKsets, nPerturbs);
    all_slope_norms = zeros(nKsets, nPerturbs);
    
    for y = 1:nPerturbs
        
        cond_index = Experimental_Data.conditions(...
            Experimental_Data.conditions == screening_conditions(y));
        rxn_indices = Experimental_Data.flux_rxns(ismember(...
            Experimental_Data.flux_rxns,uptake_rxns));
        
        uptake_values = Experimental_Data.ref_fluxes(rxn_indices,cond_index);
               
        [all_solutions{y}, all_complete_times(:,y), ...
            all_ode_warn_flags(:,y), all_slope_norms(:,y)] = ... 
            perturb_Ksets(Network_Data, All_K, All_fractions, ...
            t_interval, uptake_rxns, uptake_values);
    end
    
    for x = 1:nKsets
        
        for y = 1:nPerturbs
            
            solutions(:,y) = all_solutions{y}(:,x);
            
            fvals(x,y) = calculate_Kset_fitness(Experimental_Data, ...
                solutions(:,y), screening_conditions(y), Network_Data.S);
        
        end
        
    end
    
end