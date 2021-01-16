function [All_K,All_Rref, All_fractions, All_v_ik, count] = ...
    initial_ensemble_generation(Network_Data, nKsets, ...
    check_local_stability, uptake_rxns)
    
fprintf(1,'\n');
fprintf(1,'Generating initial ensemble of parameter sets...');
fprintf(1,'\n');

%%%%%%%%%%%%%%%%%%%%%% Define Model Characteristics %%%%%%%%%%%%%%%%%%%%%%%
v_indices = Network_Data.rxn_indices;
n_rxn_steps = ceil(((v_indices(:,2) - v_indices(:,1))+ 1)./2);
n_metabs = size(Network_Data.S,1);
nK = size(Network_Data.S_f_b,2);  
n_Rref = sum(n_rxn_steps(Network_Data.rxn_type == 1));
n_fractions = length(Network_Data.metabs_and_enzyme_complexes)- n_metabs;

%%%%%%%%%%%%%%%%%%%%%%%%%% Preallocate Variables %%%%%%%%%%%%%%%%%%%%%%%%%%
All_K = zeros(nK,nKsets);
All_Rref = zeros(n_Rref,nKsets);
All_fractions = zeros(n_fractions,nKsets);
All_v_ik = zeros(nK, nKsets);
    
n_Kset = 0;
count = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%% Sample Parameter Sets %%%%%%%%%%%%%%%%%%%%%%%%%%
    
while n_Kset < nKsets                                                       % while loop needed when screening for locally stable parameter sets to continue sampling until nKsets stables sets are sampled
        
    count = count + 1;                                                      % count informs how many parameter sets were sampled to reach nKsets; equals nKset when not checking for local stability
    
    [K, Rref, fractions, v_ik] = generate_single_set_params(Network_Data);  % sample one set of parameters (elementary kinetic parameters, reversibilities, enzyme fractions, and elementary reaction fluxes)           
        
    if check_local_stability                                                
        stable_or_nah = is_Kset_locally_stable(Network_Data, K, ...         % check if parameter set is locally stable
            fractions, uptake_rxns);      
            
        if stable_or_nah
            %%% Add parameter set to ensemble
            n_Kset = n_Kset+1;
            All_K(:,n_Kset)=K;
            All_Rref(:,n_Kset) = Rref;
            All_fractions(:,n_Kset) = fractions;
            All_v_ik(:,n_Kset) = v_ik;
        end
     
    else
        %%% Add parameter set to ensemble
        n_Kset = n_Kset+1;
        All_K(:,n_Kset)=K;
        All_Rref(:,n_Kset) = Rref;
        All_fractions(:,n_Kset) = fractions;
        All_v_ik(:,n_Kset) = v_ik;           
    end            
end

fprintf(1,'\n');
fprintf(1,'%d parameter sets sampled', count);
fprintf(1,'\n');
fprintf(1,'%d stable parameter sets generated', nKsets);
fprintf(1,'\n');
fprintf(1,'\n');
    
end
