function [Rref] = sample_reversibilities(Network_Data)

%%%%%%%%%%%%%%%%%%%%%% Define Model Characteristics %%%%%%%%%%%%%%%%%%%%%%%
v_indices = Network_Data.rxn_indices;
n_rxn_steps = ceil(((v_indices(:,2) - v_indices(:,1))+ 1)./2);
n_rev_rxns = sum(Network_Data.rxn_type == 1);
rev_rxns = find(Network_Data.rxn_type == 1);
FreeE_low = Network_Data.FreeE_range(:,1);
FreeE_high = Network_Data.FreeE_range(:,2);
WT_solution = Network_Data.WT_solution;
Rref=zeros(sum(n_rxn_steps(rev_rxns)),1);
rcount=1;

%%%%%%%%%%%%%%%%%%%%%%%% Sample reversibitilies, R %%%%%%%%%%%%%%%%%%%%%%%%

for i=1:n_rev_rxns
  
    marker = 0;
    nsteps_= n_rxn_steps(rev_rxns(i));  

    while marker == 0
        x = rand(nsteps_,1);                                                % Sample reversibility for each step in net reaction
        y = log(x);               
        
        %%% Check to see if sampled reversibilities fall within bounds:
        %%% (delG/RT)_low <= sign(Vnet)*sum(ln(R_i,j)) <= (delG/RT)_high
        bound_check = sign(WT_solution(rev_rxns(i)))*sum(y);
        
        if bound_check >= FreeE_low(rev_rxns(i)) && ...                     % Select reversibilties within bounds
                bound_check <= FreeE_high(rev_rxns(i)) 
            Rref(rcount:rcount+nsteps_-1,1)=exp(y);
            rcount=rcount+nsteps_;
            marker = 1;
            
        elseif FreeE_low(rev_rxns(i)) <= 0 && ...                           % Select reversibilties for negative delta G reactions with 0 flux
                FreeE_high(rev_rxns(i)) <= 0 && ...
                WT_solution(rev_rxns(i)) == 0
            Rref(rcount:rcount+nsteps_-1,1)=exp(y);
            rcount=rcount+nsteps_;
            marker = 1;
            
        end
    end
        
end

end