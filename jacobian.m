
function df_dx = jacobian(~,xi,k,Sr,L,T,Lo,elasticity_coeff,uptake_rxns)

    xd = T + Lo*xi;                                                         % Calculate dependent metabolite concentrations

    x = [xi; xd];                                                           % Combine independent and dependent metabolites into a single vector

    e_ij = elasticity_coeff(x,k);                                           % Calculate elasticity matrix needed to calculate jacobian of reduced model
    
    e_ij(uptake_rxns,:) = 0;                                                % uptake fluxes are constant in DAE forumlation so df_dx for uptake reactions are 0
    
    df_dx = Sr*e_ij*L;                                                      % Calculate Jacobian matrix

    df_dx = sparse(df_dx);                                                  % Convert Jacobian to sparse matrix for memory and computation time improvements

end