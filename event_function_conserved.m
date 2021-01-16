function [value,isterminal,direction] = ...
    event_function_conserved(t,xi,CPU_t,kinetic_param,S_f_b,Sr,T,Lo, ...
    uptake_rxns,uptake_values)

    Error=norm(mass_balance(0,xi,kinetic_param, ...                         % Calculate maximum deviation from dx = S*v = 0
        S_f_b,Sr,T,Lo,uptake_rxns,uptake_values));
    direction = 0;                                                          % Required for event function
    tol=1e-6;                                                               % Maximum deviation dx = S*v = 0 allowed for each metabolite before ending integration
    
    if toc(CPU_t)<300 && Error>tol                                          % If integration time is greater than 300s, end integration; 
        value=Error;
        isterminal = 0;
    elseif t <= 0.001                                                       % Avoid event not triggering if event happens on first integration step
        value=Error;
        isterminal = 0;
    else
        value=0;
        isterminal = 1;                                                     % Terminate integration if error drops below tolerance or integration time surpases 300 seconds
    end
    
end