function [K, Rref, fractions, v_ik] = generate_single_set_params(Network_Data)

%%%%%%%%%%%%%%%%%%%%%% Define Model Characteristics %%%%%%%%%%%%%%%%%%%%%%%
    n_metabs = size(Network_Data.S,1);
    S_f_b = Network_Data.S_f_b;
    A = S_f_b';
    A(A==1) = 0;
    A(A== -1) = 1;
    
%%%%%%%%%%%%%%%%%%%%%%%% Sample reversibitilies, R %%%%%%%%%%%%%%%%%%%%%%%%

    [Rref] = sample_reversibilities(Network_Data);

%%%%%%%%%%%%%%%%%%%%%% Sample enzyme fractions, Eij %%%%%%%%%%%%%%%%%%%%%%%

    fractions = sample_enzyme_fractions(Network_Data);

%%%%%%%%%%%%%%% Solve for elementary reaction rates, v_i,k %%%%%%%%%%%%%%%%
  % v_i,2j-1 = V_i,ref / (1 - R_i,j^(sign(V_i,ref)))
  % v_i,2j = V_i,ref * R_i,j^(sign(V_i,ref)) / (1 - R_i,j^(sign(V_i,ref)))

    [v_ik] = calculate_elementary_fluxes(Network_Data, Rref);

%%%%%%%%%%%%%%%%% Solve for elementary kinetic parameters %%%%%%%%%%%%%%%%%
    metab_substrate = ones(n_metabs,1);
    substrates = [metab_substrate; fractions];
    subterm = A*log(substrates);
    K = v_ik ./ exp(subterm);

end