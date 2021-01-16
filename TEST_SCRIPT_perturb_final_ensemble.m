%%% KINETIC ENSEMBLE MODELING OF E. COLI CORE METABOLISM

% Portions of code adapted from scripts provided by Liao & Maranas groups:
% 1. Rizk, M. L., and J. C. Liao. 2009. Ensemble modeling and related 
%  mathematical modeling of metabolic networks. Journal of the Taiwan 
%  Institute of Chemical Engineers 40:595-601.
% 2. Tan, Y., J. G. Rivera, C. A. Contador, J. A. Asenjo, and J. C. Liao.
%  2011. Reducing the allowable kinetic space by constructing ensemble of 
%  dynamic models with the same steady-state flux. Metab. Eng. 13:60-75.
% 3. Khodayari, A., A. R. Zomorrodi, J. C. Liao, and C. D. Maranas. 2014.
%  A kinetic model of Escherichia coli core metabolism satisfying multiple
%  sets of mutant flux data. Metab. Eng. 25:50-62.

% Portions of code adapted from previously published scripts:
%   Greene, J. L., Wäechter, A., Tyo, K. E., & Broadbelt, L. J. (2017). 
%   Acceleration Strategies to Enhance Metabolic Ensemble Modeling 
%   Performance. Biophysical Journal, 113(5), 1150-1162.

tic                                                                         % used to record calculation times                        

%%% Load Kinetic Ensemble Model of Clostridia autoethanogenum
load Final_Kinetic_Ensemble.mat
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%   SET THESE PARAMETERS PlEASE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% SET RUN PARAMETERS
                                     % INPUT uptake fluxes for CO and H2 gases.  Current values are experimental values at high biomass concentration.

  uptake_rxns = [1;2]; %[1 - CO Uptake; 2 - H2 Uptake]                        % INPUT elementary rxn index for uptake rxns constrained for DAE 
    uptake_values = [-18;-12 ];                                           % INPUT uptake fluxes for CO and H2 gases.  Current values are experimental values at high biomass concentration.

    perturbed_rxn = [1];                                                         % rxn # for enzyme to be changed; = [] if no enzyme change
    expression_level = [1];                                                      % = 1 if no enzyme change; = 0 if knockout; > 1 if overexpression; < 1 if underexpression

    t_interval = 0:5e2;                                                         % INPUT - time interval for DAE integration

    save_results = 0;                                                           % 1-save results as matlab file; 0-do not save results
    if save_results
        file_save_name = 'TEST_results';                                        % INPUT desired name of results file
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%   THANK YOU, TOUCH NOTHING ELSE, PLEASE   %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% BEGIN CODE TO SCREEN ENSEMBLE USING TRADITIONAL SCREEN %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    [solutions, complete_times, ode_warn_flags, slope_norms] = ... 
        perturb_Ksets(model, All_K_final, All_fractions_final, ...
        t_interval, uptake_rxns, uptake_values, perturbed_rxn, ...
        expression_level);
           
    time_elapsed_perturb_ensemble = toc;                                    % Prints time elapsed to screen parameter sets to command window

    %%% Generate results structure to combine ensemble predictions for a
    %%% specified test scenario
    Perturbation_Results.perturbed_rxn = perturbed_rxn;
    Perturbation_Results.expression_level = expression_level;
    Perturbation_Results.uptake_rxns = uptake_rxns;
    Perturbation_Results.uptake_values = uptake_values;
    Perturbation_Results.expression_level;
    Perturbation_Results.flux_solutions = solutions;
    Perturbation_Results.complete_times = complete_times;
    Perturbation_Results.slope_norms = slope_norms;
    Perturbation_Results.ode_warn_flags = ode_warn_flags;
    Perturbation_Results.t_interval = t_interval;
    
    clear complete_times expression_level ode_warn_flags perturbed_rxn ...
        slope_norms solutions t_interval uptake_rxns uptake_values ans
    
    if save_results
        clear save_results
        save(file_save_name)
    end

