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

rng('default')                                                              % initializes random no. generator so results are generated consistently              

%%% Indicate Elementary Model
load model_core.mat

%%% Select Reference State for Parameter Sampling
%%% High BC:    ref_condition_pub_HighBC1
%%% Med BC:     ref_condition_pub_MedBC
%%% Low BC:     ref_condition_pub_LowBC2
% cd 'Reference Conditions'\
load ref_condition_pub_HighBC1.mat
% cd ..
model.WT_solution = solution.x;
clear solution

%%% Load Experimental Data needed to calculate fitness score during
%%% screening of paramter sets
%%% Experimental_Data struct contains:
%%%     flux_rxns:      net model rxns matching experimental fluxes
%%%     conditions:     screening conditions 
%%%                     [1 - LowBC2; 2 - MedBC; 3 - HighBC1];
%%%     ref_fluxes:     experimental values for each rxn in (flux_rxns)
%%%                     for each condition in (conditions)
%%%     weighted_rxn:   rxn to weight prediction difference from
%%%                     experimental observation by (CO uptake for this
%%%                     study)
%%%     condition_names: name of each screening condition 
%%%                     (not required just for reference)    
load Experimental_Data.mat
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%   SET THESE PARAMETERS PlEASE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% SET RUN PARAMETERS
generate_Ksets = 1;                                                         % 1-generate new parameter sets; 0-use if wanting to screen previously generated parameter sets                
nKsets = 10;                                                                % INPUT number of parameter sets to generate
check_local_stability = 1;                                                  % 1-keep only sampled parameter sets that are locally stable; 0-do not screen for local stability

perturb_initial_Ksets = 1;                                                  % 1-perturb ensemble of parameter sets against condition specified below; 0-do not perturb parameter sets  
perturbation_conditions = [1 2]; %[1 - LowBC2; 2 - MedBC; 3 - HighBC1]      % INPUT conditions to test         

screen_ensemble = 1;                                                        % 1-screen ensemble against specified perturbation_conditions by ability to solve for new perturbed state
slope_norm_threshold = 1e-5;                                                % threshold of slope norm (norm of dx/dt vector for each parameter set after completing integration)
                                                                            % values close to 0 imply parameter set reached mass-balanced steady state
uptake_rxns = [1;2]; %[1 - CO Uptake; 2 - H2 Uptake]                        % INPUT elementary rxn index for uptake rxns constrained for DAE 

save_results = 1;                                                           % 1-save results as matlab file; 0-do not save results
    if save_results
        file_save_name = 'TEST_results';                                    % INPUT desired name of results file
    end
    
t_interval = 0:5e4;                                                         % INPUT - time interval for DAE integration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%   THANK YOU, TOUCH NOTHING ELSE, PLEASE   %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% BEGIN CODE TO GENERATE ENSEMBLE OF WT REF MODELS %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cd 'Ensemble Modeling Scripts'\

if generate_Ksets
    %%% Generate ensemble of nKsets parameter sets using previously
    %%% selected references state
    [All_K,All_Rref, All_fractions, All_v_ik, count] = ...
        initial_ensemble_generation(model, nKsets, ...
        check_local_stability, uptake_rxns);
    time_elapsed_parameter_set_generation = toc                             % Prints time elapsed to generate parameter sets to command window 

    Initial_Ensemble.K_sets = All_K;
    Initial_Ensemble.fractions = All_fractions;
    Initial_Ensemble.reversibilities = All_Rref;
    Initial_Ensemble.elem_fluxes = All_v_ik;
    Initial_Ensemble.locally_stable = check_local_stability;
    Initial_Ensemble.sample_count = count;

end

% All_K = matrix containing sampled elementary kinetic parameters sets
% All_Rref = matrix containing sampled reversibilities
% All_fractions = matrix containing sampled enzyme fractions
% count = number of parameter sets sampled to generate n_Ksets; should
%         equal n_Ksets if not checking for local stability

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% BEGIN CODE TO SCREEN ENSEMBLE USING TRADITIONAL SCREEN %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if perturb_initial_Ksets
    
    if screen_ensemble
        n_perturbs = length(perturbation_conditions);
        Remaining_K = All_K;
        Remaining_fractions = All_fractions;
        
        all_fvals = cell(n_perturbs,1);
        all_solutions = cell(n_perturbs,1);
        all_complete_times = cell(n_perturbs,1);
        all_ode_warn_flags = cell(n_perturbs,1);
        all_slope_norms = cell(n_perturbs,1);
        passing_ksets = cell(n_perturbs,1);
        final_kset_indices = cell(n_perturbs,1);
                
        for i = 1:n_perturbs
            
            fprintf(1,'\n');
            fprintf(1,strcat('Starting perturbation condition #', ...
                num2str(perturbation_conditions(i))));
            fprintf(1,'\n');  
            
            [all_fvals{i,1}, solutions, all_complete_times{i,1},...
                all_ode_warn_flags{i,:}, all_slope_norms{i,:}] = ...
                calculate_perturbed_steady_states(model, ...
                Experimental_Data, perturbation_conditions(i), Remaining_K, ...
                Remaining_fractions, t_interval, uptake_rxns);
                       
            passing_ksets{i,1} = find(all_slope_norms{i,1} <= slope_norm_threshold);
            
            if isempty(passing_ksets{i,1})
                fprintf(1,'\n');
                fprintf(1,'0 parameter sets remain...');
                fprintf(1,'\n');  
            else
                fprintf(1,'\n');
                fprintf(1,strcat(num2str(length(passing_ksets{i,1})),' parameter sets remain...'));
                fprintf(1,'\n');  
            end
                        
            Remaining_K = Remaining_K(:,passing_ksets{i,1});
            Remaining_fractions = Remaining_fractions(:,passing_ksets{i,1});
            all_solutions{i,1} = solutions{1,1};

        end
                              
        Screening_Results.perturbation_conditions = perturbation_conditions;
        Screening_Results.slope_norm_threshold = slope_norm_threshold;
        Screening_Results.t_interval = t_interval;
        Screening_Results.passing_ksets = passing_ksets;
        Screening_Results.solutions = all_solutions;
        Screening_Results.fitness_values = all_fvals;
        Screening_Results.complete_times = all_complete_times;
        Screening_Results.slope_norms = all_slope_norms;
        Screening_Results.ode_flags = all_ode_warn_flags;
       
        Final_Ensemble.perturbation_conditions = perturbation_conditions;
        Final_Ensemble.K_sets = Remaining_K;
        Final_Ensemble.fractions = Remaining_fractions;       
        
        final_kset_indices{n_perturbs,1,1} = passing_ksets{end,1};
        for ii = 1:length(perturbation_conditions)-1
            final_kset_indices{n_perturbs-ii,1} = passing_ksets{n_perturbs-ii,1}(final_kset_indices{n_perturbs-ii+1,1});
        end
                        
        for iii = 1:length(perturbation_conditions)
            Final_Ensemble.fitness_values(:,iii) = all_fvals{iii,1}(final_kset_indices{iii,1});
            Final_Ensemble.complete_times(:,iii) = all_complete_times{iii,1}(final_kset_indices{iii,1});
            Final_Ensemble.slope_norms(:,iii) = all_slope_norms{iii,1}(final_kset_indices{iii,1});
            Final_Ensemble.ode_flags(:,iii) = all_ode_warn_flags{iii,1}(final_kset_indices{iii,1});
            Final_Ensemble.solutions{iii,1} = all_solutions{iii,1}(:,final_kset_indices{iii,1});
        end
        
        Final_Ensemble.avg_fitness_values = mean(Final_Ensemble.fitness_values,2);
        Final_Ensemble.t_interval = t_interval;
        
        clear i ii Final_complete_times final_ensemble_sets Final_fractions ...
            Final_fvals Final_K Final_ode_warn_flags Final_slope_norms ...
            Final_slope_norms perturbation_conditions passing_ksets ...
            Final_solutions solutions iii Remaining_K Remaining_fractions ...
            final_kset_indices n_perturbs
        
    else
        
        [all_fvals, all_solutions, all_complete_times, ...
            all_ode_warn_flags, all_slope_norms] = ...
            calculate_perturbed_steady_states(model, ...
            Experimental_Data, perturbation_conditions, All_K, ...
            All_fractions, t_interval, uptake_rxns);
        
        Final_Results.perturbation_conditions = perturbation_conditions;
        Final_Results.K_sets = All_K;
        Final_Results.fractions = All_fractions;       
        Final_Results.fitness_values = all_fvals;
        Final_Results.avg_fitness_values = mean(all_fvals,2);
        Final_Results.complete_times = all_complete_times;
        Final_Results.slope_norms = all_slope_norms;
        Final_Results.ode_flags = all_ode_warn_flags;
        Final_Results.solutions = all_solutions;
        
    end
    
    time_elapsed_screen_initial_ensemble = toc;                              % Prints time elapsed to screen parameter sets to command window
end

cd ..

clear check_local_stability generate_Ksets nKsets perturb_initial_Ksets ...
    screen_ensemble uptake_rxns All_fractions All_K All_Rref All_v_ik ...
    count slope_norm_threshold all_complete_times all_fvals ...
    all_ode_warn_flags all_slope_norms all_solutions ...
    slope_norm_threshold t_interval perturbation_conditions 

if save_results
    clear save_results
    save(file_save_name)
end

% all_fvals = vector of fitness scores for each parameter set in the ensemble
% all_solutions = cell containing flux distributions for each remaining 
%   parameter set for all knockouts defined in post_screen_knockouts
% all_complete_times = cell containing computation time required to perturb
%   each remaining knockout against each specified knockout defined in
%   screening_knockouts
% all_ode_warn_flags = vector of flags for each DAE integration
%   0 - integration completed but not across the entire t_interval
%   1 - integration completed across the entire t_interval specified
%   2 - integration failed (parameter set variables set to dummy values)
%   3 - integration failed (non-real solution) and flux solution set to 0
% all_slope_norms = norm of dx/dt vector for each parameter set after
%   completing integration
