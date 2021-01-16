function [Vnet, complete_time, ode_warn_flag, norm_slope, x, t] = ...
    calculate_rates(model, t_interval, K, initial_conc,uptake_rxns, ...
    uptake_values)
    
%%%%%%%%%%%%%%%%%%%%%% Define Model Characteristics %%%%%%%%%%%%%%%%%%%%%%%
    n_rxns = size(model.S,2);
    S_f_b_conserved = model.conserved_model_info.S_f_b_conserved;
    S_f_b = model.S_f_b;
    Sr = model.conserved_model_info.Sr;
    Lo = model.conserved_model_info.Lo;
    L = model.conserved_model_info.L;
    metab_index_old_to_new = ...
        model.conserved_model_info.metab_index_old_to_new;
    metab_index_new_to_old = ...
        model.conserved_model_info.metab_index_new_to_old;
    n_independent_metabs = size(Sr,1);
    rxn_indices = model.rxn_indices;
    n_net_rxns = size(model.S,2);
    rxn_type = model.rxn_type;
       
    A = sparse(S_f_b);
    A(A == 1) = 0;
    
    %%% reorder initial concentration vector to match conservation analysis
    %%% order (independent metabolites and then dependent metabolites)
    new_initial_conc = initial_conc;
    new_initial_conc(metab_index_old_to_new) = initial_conc;
    
    %%% calculate conserved relationship totals based on initial
    %%% concentrations    
    T = new_initial_conc((n_independent_metabs+1):end) - ...
        Lo*new_initial_conc(1:n_independent_metabs);
    
    %%% model_name field is used to call model-specific elasticity
    %%% coefficient function
    model_name = model.model_name;
    elasticity_coeff = str2func(strcat('elasticity_coeff_',model_name));
    
    %%% preassign ode_warn_flag and norm_slope values; if these values are
    %%% not-reassigned, then the parameter set integration was not
    %%% successful
    ode_warn_flag = 0;
    norm_slope = 1000;
        
    try
        %%% track how long integration is running
        a = tic;        
        
        %%% apply Jacobian and Event functions for ODE solver 
        options=odeset('Jacobian',@(t,x)jacobian(t,x,K,Sr,L,T,Lo,...
            elasticity_coeff, uptake_rxns),...
            'Events',@(t,x)event_function_conserved(t,x,a,K,...
            S_f_b_conserved,Sr,T,Lo, uptake_rxns,uptake_values));
        
        %%% use stiff ode solver to calculate change in metabolite
        %%% concentration over time by integrating dx/dt equations
        [t, xi] = ode15s(@(t,x)mass_balance(t,x,K,S_f_b_conserved, ...
            Sr,T,Lo,uptake_rxns,uptake_values),t_interval,...
            new_initial_conc(1:n_independent_metabs),options);
                       
        complete_time = toc(a);
        
        %%% only calculate final integration results if solution is valid
        %%% real numbers
        if isreal(xi)          
            
            %%% note if integration completed and reached SS prior to
            %%% integrating across the entire specified time frame              
            if length(t) ~= length(t_interval)
                ode_warn_flag = 1;
            end
            
            %%% calculate net reaction fluxes from indepdent metabolite
            %%% concentrations        
            xd = repmat(T,1,size(Lo*xi',2)) + Lo*xi';                 
            x = [xi xd'];       
            x(:,metab_index_new_to_old) = x; 
            
            norm_slope = calc_norm_slope(x(end,:),K,S_f_b,uptake_rxns,...
                uptake_values);
    
            ys = repmat(x(end,:)',1,length(K));
            psub = prod(ys.^abs(A));
            vuni = diag(K*psub);
            vuni(uptake_rxns) = uptake_values;
            
            Vnet = zeros(n_net_rxns,1);
            
            for z = 1:n_net_rxns
                if rxn_type(z) == 1
                    Vnet(z) = vuni(rxn_indices(z,1)) - ...
                        vuni(rxn_indices(z,1)+1);
                    
                else
                    Vnet(z) = vuni(rxn_indices(z,1));
                    
                end
            end
       
           
        else
            Vnet = zeros(n_rxns,1);
            ode_warn_flag = 3;
        end
        
               
    catch 
        complete_time = toc(a);
        fprintf('\nError in integration. Assigning Vnet to 1000 for all rxns\n');
        Vnet = 1000*ones(n_rxns,1);
        ode_warn_flag = 2;
                

    end   
    
    
end
