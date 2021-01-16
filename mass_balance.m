function dxi = mass_balance(~,xi,k,S_f_b_conserved,Sr,T,Lo,uptake_rxns,uptake_values)

%%% calculate dependent metabolite concentrations from independent
%%% concentrations

xd = T + Lo*xi;

%%% combine independent and dependent metabolites into a single vector
x = [xi; xd]; 

%%% calculate elementary reaction fluxes
v = k;                                                      
elem_rxns = 1:length(v);
elem_rxns(uptake_rxns) = [];

for j = 1:length(elem_rxns)

    for i=1:length(x) 
        if S_f_b_conserved(i,elem_rxns(j)) < 0
            v(elem_rxns(j))=v(elem_rxns(j))* ...
                (x(i)^abs(S_f_b_conserved(i,elem_rxns(j)))); 
        end
    end
    
end

%%% DAE implementation: set uptake reaction fluxes to experimentally
%%% observed values
v(uptake_rxns) = uptake_values;


%%%%%%%%%%%%%%%%%%
dxi = Sr*v;