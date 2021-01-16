function stable_or_nah = is_Kset_locally_stable(Network_Data, K, fractions, uptake_rxns)

%%%%%%%%%%%%%%%%%%%%%% Define Model Characteristics %%%%%%%%%%%%%%%%%%%%%%%

Sr = Network_Data.conserved_model_info.Sr;
Lo = Network_Data.conserved_model_info.Lo;
metab_index_old_to_new = ...
    Network_Data.conserved_model_info.metab_index_new_to_old;
no_independent_metabs = size(Sr,1);
L = Network_Data.conserved_model_info.L;
no_metabs = size(Network_Data.S,1);       

%%% reorder concentration vector to match conserved information metabolite
%%% order
initial_conc = [ones(no_metabs,1); fractions];
new_initial_conc = initial_conc;
new_initial_conc(metab_index_old_to_new) = initial_conc;

%%% calculate conserved relationship totals based on initial concentrations
T = new_initial_conc((no_independent_metabs+1):end) - ...
    Lo*new_initial_conc(1:no_independent_metabs);

elasticity_coeff = str2func(strcat('elasticity_coeff_', ...
    Network_Data.model_name));

%%%%%%%%%%%%%% Calculate Jacobian with sampled parameter set %%%%%%%%%%%%%%
%%%%%%%%%%%%%%% @ reference steady state flux distribution %%%%%%%%%%%%%%%%

df_dx = full(jacobian(0,new_initial_conc(1:no_independent_metabs), ...
    K,Sr,L,T,Lo, elasticity_coeff, uptake_rxns));
eigen_values = eig(df_dx);                                                  % Calculate eigenvalues of Jacobian
real_parts = real(eigen_values);                                            % Identify real part values of eigenvalues
    
    if max(real_parts) <= -0.0000001                                        % If all real parts are negative, the parameter set is marked as locally stable
        stable_or_nah = 1;
    else
        stable_or_nah = 0;
    
    end

end
