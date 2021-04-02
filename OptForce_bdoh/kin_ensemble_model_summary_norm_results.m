
optForceSets = load('./Final/Backup/optForceSets.mat');
typeRegOptForceSets = load('./Final/Backup/typeRegOptForce.mat'); 

%selection of first order OptForce rxns that are in core model
load('./Edited_Results/First_Order_Core.mat')
%selection of second order OptForce rxns that are in core model
load('./Final/Second_Order_Core.mat')
%selection of third order OptForce rxns that are in core model
load('./Final/Third_Order_Core.mat')

model = readCbModel('./Final/model_core.mat');
for i=1:3
    
    First_order_results(:,27,i)= Base_BDOH';
    Second_order_results(:,27,i)=Base_BDOH';
    Third_order_results(:,27,i)=Base_BDOH';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting section for kinetic ensembel model results n first order
%interventions suggested by Optforce
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('./Edited_Results/Normalized_first.mat')

load('./Edited_Results/First_Order_Results.mat')
First_order_Rxn = model.rxns(First_Order_Core(:,1));
First_order_Intervent = First_Order_Core(:,2)';

%clear 
clf

%create a panel of plots distributed in two rows and three columns. The third
%number specifies the position that the plot is taking in the panel.

%Plot the BDO reaction flux according to each intervention and expression
%regulation level 1 
subplot(2,3,1)
h1 = heatmap(Normalized_Results_First(:,:,1));
h1.XDisplayLabels = First_order_Rxn;
h1ip = get(h1,'InnerPosition');

%Plot the type of intevention of each reaction suggested by Optforce
subplot(2,3,4)
h2 = heatmap(First_order_Intervent);
h2.XDisplayLabels = First_order_Rxn;
set(h2, 'InnerPosition', [h1ip(1) 0.45 h1ip(3)-0.055 0.025]); 

%Plot the BDO reaction flux according to each intervention and expression
%regulation level 2
subplot(2,3,2)
h3 = heatmap(Normalized_Results_First(:,:,2));
h3.XDisplayLabels = First_order_Rxn;
h3ip = get(h3,'InnerPosition');

subplot(2,3,5)
h4 = heatmap(First_order_Intervent);
h4.XDisplayLabels = First_order_Rxn;
set(h4, 'InnerPosition', [h3ip(1) 0.45 h3ip(3)-0.055 0.025]); 

%Plot the BDO reaction flux according to each intervention and expression
%regulation level 3
subplot(2,3,3)
h5 = heatmap(Normalized_Results_First(:,:,3));
h5.XDisplayLabels = First_order_Rxn;
h5ip = get(h5,'InnerPosition');

subplot(2,3,6)
h6 = heatmap(First_order_Intervent);
h6.XDisplayLabels = First_order_Rxn;
set(h6, 'InnerPosition', [h5ip(1) 0.45 h5ip(3)-0.055 0.025]); 

%save the figure as maually to preserve the format. Otherwise
%saveas(gcf,'kem.norm.first.order.intervenions.jpeg');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting section for kinetic ensembel model results on second order
%interventions suggested by Optforce
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('./Edited_Results/Normalized_Second.mat')

Second_order_Rxn_A = model.rxns(Second_Order_Core(:,1));
Second_order_Rxn_B = model.rxns(Second_Order_Core(:,3));
Second_order_Intervent_A = Second_Order_Core(:,2)';
Second_order_Intervent_B = Second_Order_Core(:,4)';

%clear 
clf

%create a panel of plots distributed in two rows and three columns. The third
%number specifies the position that the plot is taking in the panel.

%Plot the BDO reaction flux according to each intervention and expression
%regulation level 1 
subplot(3,3,1)
h1 = heatmap(Normalized_Results_Second(:,:,1));
h1ip = get(h1,'InnerPosition');

%Plot the type of intevention of each reaction A suggested by Optforce
subplot(3,3,4)
h2 = heatmap(Second_order_Intervent_A);
h2.ColorLimits = [0 2];
h2.XDisplayLabels = Second_order_Rxn_A;
set(h2, 'InnerPosition', [h1ip(1) 0.65 h1ip(3)-0.055 0.025]); 

%Plot the type of intevention of each reaction B suggested by Optforce
subplot(3,3,7)
h3 = heatmap(Second_order_Intervent_B);
h3.ColorLimits = [0 2];
h3.XDisplayLabels = Second_order_Rxn_B;
set(h3, 'InnerPosition', [h1ip(1) 0.55 h1ip(3)-0.055 0.025]);

%Plot the BDO reaction flux according to each intervention and expression
%regulation level 2 
subplot(3,3,2)
h4 = heatmap(Normalized_Results_Second(:,:,2));
h4ip = get(h4,'InnerPosition');

%Plot the type of intevention of each reaction A suggested by Optforce
subplot(3,3,5)
h5 = heatmap(Second_order_Intervent_A);
h5.ColorLimits = [0 2];
h5.XDisplayLabels = Second_order_Rxn_A;
set(h5, 'InnerPosition', [h4ip(1) 0.65 h4ip(3)-0.055 0.025]); 

%Plot the type of intevention of each reaction B suggested by Optforce
subplot(3,3,8)
h6 = heatmap(Second_order_Intervent_B);
h6.ColorLimits = [0 2];
h6.XDisplayLabels = Second_order_Rxn_B;
set(h6, 'InnerPosition', [h4ip(1) 0.55 h4ip(3)-0.055 0.025]); 


%Plot the BDO reaction flux according to each intervention and expression
%regulation level 3 
subplot(3,3,3)
h7 = heatmap(Normalized_Results_Second(:,:,3));
h7ip = get(h7,'InnerPosition');

%Plot the type of intevention of each reaction A suggested by Optforce
subplot(3,3,6)
h8 = heatmap(Second_order_Intervent_A);
h8.ColorLimits = [0 2];
h8.XDisplayLabels = Second_order_Rxn_A;
set(h8, 'InnerPosition', [h7ip(1) 0.65 h7ip(3)-0.055 0.025]); 

%Plot the type of intevention of each reaction B suggested by Optforce
subplot(3,3,9)
h9 = heatmap(Second_order_Intervent_B);
h9.ColorLimits = [0 2];
h9.XDisplayLabels = Second_order_Rxn_B;
set(h9, 'InnerPosition', [h7ip(1) 0.55 h7ip(3)-0.055 0.025]); 

%save as figure from the automatically generated figure


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting section for kinetic ensembel model results on third order
%interventions suggested by Optforce
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('./Edited_Results/Normalized_Third.mat')

Third_order_Rxn_A = model.rxns(Third_Order_Core(:,1));
Third_order_Rxn_B = model.rxns(Third_Order_Core(:,3));
Third_order_Rxn_C = model.rxns(Third_Order_Core(:,5));

Third_order_Intervent_A = Third_Order_Core(:,2)';
Third_order_Intervent_B = Third_Order_Core(:,4)';
Third_order_Intervent_C = Third_Order_Core(:,6)';

%clear 
clf

%create a panel of plots distributed in two rows and three columns. The third
%number specifies the position that the plot is taking in the panel.

%Plot the BDO reaction flux according to each intervention and expression
%regulation level 1 
subplot(4,3,1)
h1 = heatmap(Normalized_Results_Third(:,:,1));
h1ip = get(h1,'InnerPosition');

%Plot the type of intevention of each reaction A suggested by Optforce
subplot(4,3,4)
h2 = heatmap(Third_order_Intervent_A);
h2.ColorLimits = [0 2];
h2.XDisplayLabels = Third_order_Rxn_A;
set(h2, 'InnerPosition', [h1ip(1) 0.65 h1ip(3)-0.055 0.025]); 

%Plot the type of intevention of each reaction B suggested by Optforce
subplot(4,3,7)
h3 = heatmap(Third_order_Intervent_B);
h3.ColorLimits = [0 2];
h3.XDisplayLabels = Third_order_Rxn_B;
set(h3, 'InnerPosition', [h1ip(1) 0.55 h1ip(3)-0.055 0.025]);

%Plot the type of intevention of each reaction C suggested by Optforce
subplot(4,3,10)
h4 = heatmap(Third_order_Intervent_C);
h4.ColorLimits = [0 2];
h4.XDisplayLabels = Third_order_Rxn_C;
set(h4, 'InnerPosition', [h1ip(1) 0.45 h1ip(3)-0.055 0.025]);

%Plot the BDO reaction flux according to each intervention and expression
%regulation level 2 
subplot(4,3,2)
h5 = heatmap(Normalized_Results_Third(:,:,2));
h5ip = get(h5,'InnerPosition');

%Plot the type of intevention of each reaction A suggested by Optforce
subplot(4,3,5) 
h6 = heatmap(Third_order_Intervent_A);
h6.ColorLimits = [0 2];
h6.XDisplayLabels = Third_order_Rxn_A;
set(h6, 'InnerPosition', [h5ip(1) 0.65 h5ip(3)-0.055 0.025]); 

%Plot the type of intevention of each reaction B suggested by Optforce
subplot(4,3,8)
h7 = heatmap(Third_order_Intervent_B);
h7.ColorLimits = [0 2];
h7.XDisplayLabels = Third_order_Rxn_B;
set(h7, 'InnerPosition', [h5ip(1) 0.55 h5ip(3)-0.055 0.025]);

%Plot the type of intevention of each reaction C suggested by Optforce
subplot(4,3,11)
h8 = heatmap(Third_order_Intervent_C);
h8.ColorLimits = [0 2];
h8.XDisplayLabels = Third_order_Rxn_C;
set(h8, 'InnerPosition', [h5ip(1) 0.45 h5ip(3)-0.055 0.025]);

%Plot the BDO reaction flux according to each intervention and expression
%regulation level 3
subplot(4,3,3)
h9 = heatmap(Normalized_Results_Third(:,:,3));
h9ip = get(h9,'InnerPosition');

%Plot the type of intevention of each reaction A suggested by Optforce
subplot(4,3,6)
h10 = heatmap(Third_order_Intervent_A);
h10.ColorLimits = [0 2];
h10.XDisplayLabels = Third_order_Rxn_A;
set(h10, 'InnerPosition', [h9ip(1) 0.65 h9ip(3)-0.055 0.025]); 

%Plot the type of intevention of each reaction B suggested by Optforce
subplot(4,3,9)
h11 = heatmap(Third_order_Intervent_B);
h11.ColorLimits = [0 2];
h11.XDisplayLabels = Third_order_Rxn_B;
set(h11, 'InnerPosition', [h9ip(1) 0.55 h9ip(3)-0.055 0.025]);

%Plot the type of intevention of each reaction C suggested by Optforce
subplot(4,3,10)
h12 = heatmap(Third_order_Intervent_C);
h12.ColorLimits = [0 2];
h12.XDisplayLabels = Third_order_Rxn_C;
set(h12, 'InnerPosition', [h9ip(1) 0.45 h9ip(3)-0.055 0.025]);

%save as figure from the automatically generated figure


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Statistical analysis 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
   

%% Apply statistical test for the results on first order interventions

%when the regulation of the expression of the ensyme is low (1)
[p,tbl,stats] = anova1(Normalized_Results_First(:,:,1));
[c,m,h,gnames] = multcompare(stats);

%when the regulation of the expression of the ensyme is medium (2)
[p,tbl,stats] = anova1(Normalized_Results_First(:,:,2));
[c,m,h,gnames] = multcompare(stats);

%when the regulation of the expression of the ensyme is high (3)
[p,tbl,stats] = anova1(Normalized_Results_First(:,:,3));
[c,m,h,gnames] = multcompare(stats);

%The BDO flux obtained by down-regulation of Rnf (regulatiry level 3) is
%significanty higher than the majority of the remaining interventions 


%% Apply the statistical test for the results on second order interventions

%when the regulation of the expression of the ensyme is low (1)
[p,tbl,stats] = anova1(Normalized_Results_Second(:,:,1));
[c,m,h,gnames] = multcompare(stats);

%when the regulation of the expression of the ensyme is medium (2)
[p,tbl,stats] = anova1(Normalized_Results_Second(:,:,2));
[c,m,h,gnames] = multcompare(stats);

%when the regulation of the expression of the ensyme is high (3)
[p,tbl,stats] = anova1(Normalized_Results_Second(:,:,3));
[c,m,h,gnames] = multcompare(stats);

%No intervention significantly differs from others at whatever regulatory
%level

%% Apply the statistical test for the results on third order interventions

%when the regulation of the expression of the ensyme is low (1)
[p,tbl,stats] = anova1(Normalized_Results_Third(:,:,1));
[c,m,h,gnames] = multcompare(stats);

%when the regulation of the expression of the ensyme is medium (2)
[p,tbl,stats] = anova1(Normalized_Results_Third(:,:,2));
[c,m,h,gnames] = multcompare(stats);

%when the regulation of the expression of the ensyme is high (3)
[p,tbl,stats] = anova1(Normalized_Results_Third(:,:,3));
[c,m,h,gnames] = multcompare(stats);

%No intervention significantly differs from others at whatever regulatory
%level