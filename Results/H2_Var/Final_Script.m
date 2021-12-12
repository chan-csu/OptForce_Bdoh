load('EXP_Struct.mat')

for i=1:size(EXP_Struct,2)
    Uptakes.rxns{i}=EXP_Struct(i).Name;
    Uptakes.values{i}=EXP_Struct(i).EXP1;                                     %If you want another experiment, add that to EXP_Struct, and change EXP1 accordingly in this line
    Uptakes.type{i}='b';
end

load('Metaclau.mat')

[optForceSets, posOptForceSets, typeRegOptForceSets, flux_optForceSets,Eff_out] ...
    =Opt_find(metaclau,Uptakes)

load('Mapping vect.mat')
load('model_core.mat')
load Final_Kinetic_Ensemble.mat
%% Generation of perturbation vectors for the Core/Kinetic model
Temp=optForceSets{1,1};
j=0
for i=1:size(optForceSets{1,1},1)
    [a,b]=ismember(Temp{i},Mapping_Vect(:,1))
    if a
        j=j+1
        First_Order_Core(j,1)=find(ismember(model.rxns(:),(Mapping_Vect(find(ismember(Mapping_Vect(:,1),Temp{i})),2))))
        Temp2=typeRegOptForceSets{1,1}
        if strcmp(Temp2{i},'knockout')
            First_Order_Core(j,2)=0;
        elseif strcmp(Temp2{i},'downregulation')
            First_Order_Core(j,2)=1;
        else
            First_Order_Core(j,2)=2                                         %I made a convention here: Knockout:0 Downreg:1 Upreg=1 for convenience 
        end
    end
end

Temp=optForceSets{1,2};
j=0
for i=1:size(optForceSets{1,2},1)
    [a,b]=ismember(Temp{i,1},Mapping_Vect(:,1))
    [c,d]=ismember(Temp{i,2},Mapping_Vect(:,1))
    if a && c
        j=j+1
        Second_Order_Core(j,1)=find(ismember(model.rxns,(Mapping_Vect(find(ismember(Mapping_Vect(:,1),Temp{i,1})),2))))
        Second_Order_Core(j,3)=find(ismember(model.rxns,(Mapping_Vect(find(ismember(Mapping_Vect(:,1),Temp{i,2})),2))))
        Temp2=typeRegOptForceSets{1,2}
        if strcmp(Temp2{i,1},'knockout')
            Second_Order_Core(j,2)=0;
        elseif strcmp(Temp2{i,1},'downregulation')
            Second_Order_Core(j,2)=1;
        else
            Second_Order_Core(j,2)=2                                         %I made a convention here: Knockout:0 Downreg:1 Upreg=1 for convenience 
        end
        
        if strcmp(Temp2{i,2},'knockout')
            Second_Order_Core(j,4)=0;
        elseif strcmp(Temp2{i,2},'downregulation')
            Second_Order_Core(j,4)=1;
        else
            Second_Order_Core(j,4)=2                                         
        end
    end
end

Temp=optForceSets{1,3};
j=0
for i=1:size(optForceSets{1,3},1)
    [a,b]=ismember(Temp{i,1},Mapping_Vect(:,1))
    [c,d]=ismember(Temp{i,2},Mapping_Vect(:,1))
    [e,f]=ismember(Temp{i,3},Mapping_Vect(:,1))
    if a && c && e
        j=j+1
        Third_Order_Core(j,1)=find(ismember(model.rxns,(Mapping_Vect(find(ismember(Mapping_Vect(:,1),Temp{i,1})),2))))
        Third_Order_Core(j,3)=find(ismember(model.rxns,(Mapping_Vect(find(ismember(Mapping_Vect(:,1),Temp{i,2})),2))))
        Third_Order_Core(j,5)=find(ismember(model.rxns,(Mapping_Vect(find(ismember(Mapping_Vect(:,1),Temp{i,3})),2))))
      
        Temp2=typeRegOptForceSets{1,3}
        if strcmp(Temp2{i,1},'knockout')
            Third_Order_Core(j,2)=0;
        elseif strcmp(Temp2{i,1},'downregulation')
            Third_Order_Core(j,2)=1;
        else
            Third_Order_Core(j,2)=2                                         %I made a convention here: Knockout:0 Downreg:1 Upreg=1 for convenience 
        end
        
        if strcmp(Temp2{i,2},'knockout')
            Third_Order_Core(j,4)=0;
        elseif strcmp(Temp2{i,2},'downregulation')
            Third_Order_Core(j,4)=1;
        else
            Third_Order_Core(j,4)=2                                         
        end
        
         if strcmp(Temp2{i,3},'knockout')
            Third_Order_Core(j,6)=0;
        elseif strcmp(Temp2{i,2},'downregulation')
            Third_Order_Core(j,6)=1;
        else
            Third_Order_Core(j,6)=2                                         
         end
        
         
    end
end
            
        
%% First Order Interventions
Upregulation_Vect=[3 5 10];
Downregulation_Vect=[0.3 0.2 0.1];
Knockout_Vect=[0 0 0];
uptake_rxns = [1;2]; %[1 - CO Uptake; 2 - H2 Uptake]                        % INPUT elementary rxn index for uptake rxns constrained for DAE 
uptake_values = [Uptakes.values{1};Uptakes.values{2} ];                                           % INPUT uptake fluxes for CO and H2 gases.  Current values are experimental values at high biomass concentration.

perturbed_rxn = [1];                                                         % rxn # for enzyme to be changed; = [] if no enzyme change
expression_level = [1];                                                      % = 1 if no enzyme change; = 0 if knockout; > 1 if overexpression; < 1 if underexpression
t_interval = 0:5e2;                                                         % INPUT - time interval for DAE integration

    save_results = 0;                                                           % 1-save results as matlab file; 0-do not save results
    if save_results
        file_save_name = 'TEST_results';                                        % INPUT desired name of results file
    end



    
    [solutions, complete_times, ode_warn_flags, slope_norms] = ... 
        perturb_Ksets(model, All_K_final, All_fractions_final, ...
        t_interval, uptake_rxns, uptake_values, perturbed_rxn, ...
        expression_level);
           
Base_BDOH=solutions(62,:);
T=array2table(solutions);
T.Properties.RowNames=model.rxns(:)
writetable(T,'Kinetic_Model_Fluxes.xlsx','Sheet','Base','WriteRowNames',true)
for i=1:size(First_Order_Core,1)
    for j=1:length(Upregulation_Vect)
        if First_Order_Core(i,2)==0
            perturbed_rxn=[First_Order_Core(i,1)];
            expression_level = Knockout_Vect(j);
        elseif First_Order_Core(i,2)==1
            perturbed_rxn=[First_Order_Core(i,1)];
            expression_level = Downregulation_Vect(j);
        else
            perturbed_rxn=[First_Order_Core(i,1)];
            expression_level = Upregulation_Vect(j);
        end
    [solutions, complete_times, ode_warn_flags, slope_norms] = ... 
        perturb_Ksets(model, All_K_final, All_fractions_final, ...
        t_interval, uptake_rxns, uptake_values, perturbed_rxn, ...
        expression_level);
        First_order_results(:,i,j)=solutions(62,:)
        T=array2table(solutions);
        T.Properties.RowNames=model.rxns(:)
        writetable(T,'Kinetic_Model_Fluxes.xlsx','Sheet',strcat('First_order',num2str(i)),'WriteRowNames',true)

    end
end
for i=1:size(Second_Order_Core,1)
    for j=1:length(Upregulation_Vect)
        if Second_Order_Core(i,2)==0
            perturbed_rxn(1)=[Second_Order_Core(i,1)];
            expression_level(1) = Knockout_Vect(j);
        elseif Second_Order_Core(i,2)==1
            perturbed_rxn(1)=[Second_Order_Core(i,1)];
            expression_level(1) = Downregulation_Vect(j);
        else
            perturbed_rxn(1)=[Second_Order_Core(i,1)];
            expression_level(1)= Upregulation_Vect(j)
        end
        
        if Second_Order_Core(i,4)==0
            perturbed_rxn(2)=[Second_Order_Core(i,3)];
            expression_level(2) = Knockout_Vect(j);
        elseif Second_Order_Core(i,2)==1
            perturbed_rxn(2)=[Second_Order_Core(i,3)];
            expression_level(2) = Downregulation_Vect(j);
        else
            perturbed_rxn(2)=[Second_Order_Core(i,3)];
            expression_level(2) = Upregulation_Vect(j);
        end
    [solutions, complete_times, ode_warn_flags, slope_norms] = ... 
        perturb_Ksets(model, All_K_final, All_fractions_final, ...
        t_interval, uptake_rxns, uptake_values, perturbed_rxn, ...
        expression_level);
        Second_order_results(:,i,j)=solutions(62,:)
        T=array2table(solutions);
        T.Properties.RowNames=model.rxns(:)
        writetable(T,'Kinetic_Model_Fluxes.xlsx','Sheet',strcat('Second_order',num2str(i)),'WriteRowNames',true)

    end
end


for i=1:size(Third_Order_Core,1)
    for j=1:length(Upregulation_Vect)
        if Third_Order_Core(i,2)==0
            perturbed_rxn(1)=[Third_Order_Core(i,1)];
            expression_level(1) = Knockout_Vect(j)
        elseif Third_Order_Core(i,2)==1
            perturbed_rxn(1)=[Third_Order_Core(i,1)];
            expression_level(1) = Downregulation_Vect(j)
        else
            perturbed_rxn(1)=[Third_Order_Core(i,1)];
            expression_level(1)= Upregulation_Vect(j)
        end
        
        if Third_Order_Core(i,4)==0
            perturbed_rxn(2)=[Third_Order_Core(i,3)];
            expression_level(2) = Knockout_Vect(j)
        elseif Third_Order_Core(i,4)==1
            perturbed_rxn(2)=[Third_Order_Core(i,3)];
            expression_level(2) = Downregulation_Vect(j)
        else
            perturbed_rxn(2)=[Third_Order_Core(i,3)];
            expression_level(2) = Upregulation_Vect(j)
        end
        
        if Third_Order_Core(i,6)==0
            perturbed_rxn(3)=[Third_Order_Core(i,5)];
            expression_level(3) = Knockout_Vect(j)
        elseif Third_Order_Core(i,6)==1
            perturbed_rxn(3)=[Third_Order_Core(i,5)];
            expression_level(3) = Downregulation_Vect(j)
        else
            perturbed_rxn(3)=[Third_Order_Core(i,5)];
            expression_level(3) = Upregulation_Vect(j)
        end
        
        
    [solutions, complete_times, ode_warn_flags, slope_norms] = ... 
        perturb_Ksets(model, All_K_final, All_fractions_final, ...
        t_interval, uptake_rxns, uptake_values, perturbed_rxn, ...
        expression_level);
        Third_order_results(:,i,j)=solutions(62,:)
        T=array2table(solutions);
        T.Properties.RowNames=model.rxns(:)
        writetable(T,'Kinetic_Model_Fluxes.xlsx','Sheet',strcat('Third_order',num2str(i)),'WriteRowNames',true)


    end
end
     

for i=1:18
Normalized_Results_First(i,:,:)=First_order_results(i,:,:)./Base_BDOH(i)
end
for i=1:18
Normalized_Results_Second(i,:,:)=Second_order_results(i,:,:)./Base_BDOH(i)
end
for i=1:18
Normalized_Results_Third(i,:,:)=Third_order_results(i,:,:)./Base_BDOH(i)
end
