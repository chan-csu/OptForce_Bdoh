function [v_ik] = calculate_elementary_fluxes(Network_Data, Rref)
     
%%%%%%%%%%%%%%%%%%%%%% Define Model Characteristics %%%%%%%%%%%%%%%%%%%%%%%

rxn_type = Network_Data.rxn_type;  
v_indices = Network_Data.rxn_indices;
n_rxn_steps = ceil(((v_indices(:,2) - v_indices(:,1))+ 1)./2);
n_rev_rxns = sum(rxn_type == 1);
rev_rxns = find(rxn_type == 1);
n_elem_rxns = length(Network_Data.rxns_f_b);
n_transport_rxns = sum(rxn_type == 2);
transport_rxns = find(rxn_type == 2);
biomass_rxns = find(rxn_type == 4);
n_biomass_rxns = sum(rxn_type == 4);

rVnet = Network_Data.WT_solution;
rVnet(rVnet > -1e-10 & rVnet < 1e-10) = 0;                                  % set small rxn flux values to 0 for v_ik calculations

v_ik = zeros(n_elem_rxns,1);                                                % preallocate elementary flux vector

rcount = 0;

%%%%%%%%%%%%%%% Solve for elementary reaction rates, v_i,k %%%%%%%%%%%%%%%%
  % v_i,2j-1 = V_i,ref / (1 - R_i,j^(sign(V_i,ref)))
  % v_i,2j = V_i,ref * R_i,j^(sign(V_i,ref)) / (1 - R_i,j^(sign(V_i,ref)))

for i = 1:n_rev_rxns                                                        % calculate elementary fluxes for one reaction at a time
    
    n_steps = n_rxn_steps(rev_rxns(i));
    v_index = v_indices(rev_rxns(i),1):1:v_indices(rev_rxns(i),2);
    sign_rVnet = sign(rVnet(rev_rxns(i)));
    Vnet_i = rVnet(rev_rxns(i));
    
    if Vnet_i ~= 0                                                        
        
        for j = 1:n_steps

            R_step = Rref(rcount + j);
            v_ik(v_index(2*j-1)) = Vnet_i/(1 - R_step^(sign_rVnet));        % calculate forward elementary fluxes from WT rxn flux distribution and sampled reaction revseribilities
            v_ik(v_index(2*j)) = ...                                        % calculate reverse elementary fluxes from WT rxn flux distribution and sampled reaction revseribilities
                Vnet_i*R_step^(sign_rVnet)/(1 - R_step^(sign_rVnet));
            
            if v_ik(v_index(2*j-1)) < 0                                     % if reaction is reversed, switch forward and backward reaction rates
                switch_forward = v_ik(v_index(2*j-1));
                switch_backward = v_ik(v_index(2*j));
                v_ik(v_index(2*j-1)) = -switch_backward;
                v_ik(v_index(2*j))= -switch_forward;
            end  
            
        end
        
    else                                                                    % if net reaction flux is 0, randomly sample elementary reaction fluxes within range of net reaction flux distribution
        
        for j = 1:n_steps
            
            v_ik(v_index(2*j-1)) = rand(1,1)*abs(max(rVnet));               % sample forward elementary reaction flux
            v_ik(v_index(2*j)) = v_ik(v_index(2*j-1));                      % set reverse elementary reaction flux to forward value so net is equal to 0
            
        end       
           
    end
    rcount=rcount+n_steps;
end
    
for m = 1:n_transport_rxns                                                  % set elementary fluxes for transport reactions to net reaction flux value
    
    v_ik(v_indices(transport_rxns(m),1)) = rVnet(transport_rxns(m));
    
end

for b = 1:n_biomass_rxns
    v_ik(v_indices(biomass_rxns(b),1):v_indices(biomass_rxns(b),2)) = ...
        rVnet(biomass_rxns(b)); 
    
end

end
