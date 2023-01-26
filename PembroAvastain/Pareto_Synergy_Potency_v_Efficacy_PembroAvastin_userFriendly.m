%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
%% Multi-objective synergy optimization by Jana Gevertz (9/7/2022)     %%
% - The code first establishes the pembro+beva model, its parameters,   %
%   initial conditions, and dosing protocol.                            %
% - With all model parameters except the dose of drug fixed, we sweep   %
%   over a range of monotherapy doses (d1 and d2 values)                %
% - The code calculates the tumor growth inhibition rate RELATIVE to    %
%   control case at a fixed time point tf using the formula:            %
%    TGI = (V(control)-V(drug))/V(control)                              %
% - Multiple combination indices CI are computed, each using a          %
%   different definition of additivity (no interaction):                %
%   1) Synergy of efficacy using highest single agent (HSA):            %
%       CI(HSA) = max[E(d1),E(d2)]/E(d1,d2)                             %
%   2) Synergy of efficacy using Bliss                                  %
%       CI(Bliss) = [E(d1)+E(d2)-E(d1)E(d2)]/E(d1,d2)                   %
%   3) Synergy of potency using novel method lowest single dose (LSD)   %
%       CI(LSD) = (d1+d2)/min(D1,D2)                                    %
%      where E(D1,0)=E(0,D2)=E(d1,d2). I currently measure the potency  %
%      of combination as sum of drug doses, which (probably) assumes    %
%      the two drugs equally contribute to toxicity. If one drug is     %
%      more toxic, this input could be changed to an appropriately-     %
%      weighted sum. Note: while this is still implemented, we will     %
%      instead using a relative LSD score:                              %
%       CI(LSDrel) = ((d1/EC50_1)+(d2/EC50_2))/min(D1/EC50_1,D2/EC50_2) %
%      where EC50_i refers to the DOSE of drug i that results in 50%    %
%      inhibition relative to control case                              %
%   4) Synergy of potency using Loewe's additivity formula:             %
%       CI(Loewes) = (d1/D1)+(d2/D2)                                    %
%      where CI = 1 is additive, CI<1 is synergy, CI>1 is antagonism    %
%      (Note: could use generalized Loewe's CI instead: see code        %
%      Understanding_Loewes/Loewes_GeneralizedLoewes_CI.m               %
% - The code plots:                                                     %         
%   Figure 1: heatmap of drug 1 monotherapy TGI relative to control     %          
%   Figure 2: heatmap of drug 2 monotherapy TGI relative to control     % 
%   Figure 3: heatmap of combination TGI relative to control            %
%   Figure 4: heatmap of synergy of efficacy defined by HSA             %
%   Figure 5: heatmap of synergy of efficacy defined by Bliss           %
%   Figure 6: heatmap of cumulative relative dose (relative to EC50)    %
%   Figure 7: heatmap of synergy of potency defined by relative LSD     %
%   Figure 8: heatmap of synergy of potency defined by Loewe           %
%   Figure 9: Contour plot of isobolograms (curves of constant         %
%              combination efficacy)                                    %
%   Figure 10: Pareto front rel. LSD-HSA multi-synergy criterion        %
%   Figure 11: Doses that fall on relative LSD-HSA Pareto front,        %
%              displayed on TGI heatmap                                 %
%   Figure 12: Pareto front rel. LSD-Bliss multi-synergy criterion space%
%   Figure 13: Doses that fall on relative LSD-Bliss Pareto front       %
%              displayed on TGI heatmap                                 %
%   Figure 14: Pareto front Loewe-HSA multi-synergy criterion space     %
%   Figure 15: Doses that fall on Loewe-HSA Pareto front                %
%              displayed on TGI heatmap                                 %
%   Figure 16: Pareto front Loewe-Bliss multi-synergy criterion space   %
%   Figure 17: Doses that fall on Loewe-Bliss Pareto front              %
%              displayed on TGI heatmap                                 %
%   Figure 18: (Added 10/5) On one plot, shows all doses that appear on %
%              all Pareto fronts, with the number shown indicating how  %
%              many Pareto fronts the dose fell on                      %
% - All output is written to PointsOnPareto, including a list of all    %
%   doses on each Pareto front, and a list of all doses that appear on  %
%   multiple Pareto fronts                                              %
% - Run time (num_pts = 30, cell line = 2): 7450.287158 seconds         %
% - Updates on 10/5:                                                    %
%   1) Include experimental combination dose on heatmap that shows the  %
%      Pareto optimal doses. For H1299, had to expand max monotherapy   %
%      range to EC66 from EC60, or else experimental dose was not in    %
%      dose range considered.                                           %
%   2) Introduce variables d1_spacing and d2_spacing so can repeat      %
%      analysis for other protocols                                     %
%% - The code is readily adaptable to other models. User would need to %%
%   appropriately modify set_parameters_IC_protocol() to define the     %
%   model parameters, initial conditions, and dosing protocol. User     %
%   would also need to determine the range of doses to consider for     %
%   thier model. Finally, the user needs to ensure their code that      %
%   solves the DE just returns the tumor size at the final time point.  %
%   If all these elements are in place, all the combination indices     %
%   will automatically be computed, and the plots will all be made.     %
%   As of 10/5 update, also can modify spacing of doses.                %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; clc; close all; tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Declare model parameters, ICs and dosing protocol                   %%
%% **User can re-define this function to correspond to any model**     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d1_spacing = 3; % frequency of drug 1 (experiments: given every 3 days)
d2_spacing = 3; % frequency of drug 2 (experiments: given every 3 days)
num_pts = 30; % number of discrete doses to consider over domain of each dose
[p,in_cond,tmax,tStart1,tStart2,d1,d2,d1_search,d2_search,max_d1,min_d2,max_d2] = ...
    set_parameters_ICs_protocol(d1_spacing,d2_spacing,num_pts); 
if p.cell_line == 1
    path = 'Pembro_Bevacizumab/Output_H1299_ChangeSearch';
    fprintf('Cell line %d = H1299\n',p.cell_line);
elseif p.cell_line == 2
    path = 'Pembro_Bevacizumab/Output_A549_ChangeSearch';
    fprintf('Cell line %d = A549\n',p.cell_line);
end
if exist(path, 'dir') ~= 7
    mkdir(path)
end
sdiary = [path '/PointsOnPareto'];
diary(sdiary); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Range of doses to consider is now set for each drug (for each cell  %%
%% line. This does need to be manually determined for a given model.   %%
%% Here the smallest dose to consider will be the EC5, so that         %%
%% determines min_d1, min_d2. The largest monotherapy dose to consider %%
%% is constrained by the largest combination efficacy. Here went with  %%
%% EC60 for largest monotherapy dose, so that determines max_d1,       %%
%% max_d2. For each cell line, the max efficacy of the combination     %%
%% therapy is 94-95%. On 10/5, changes to EC66 as max for H1299 so     %%
%% experimental dose was in dose range searched. Max efficacy of combo %%
%% remains in 94-95% range.                                            %%
%% **User must determine appropriate dose range for their model**      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d1_experiment = 10; % actual dose of drug 1
d2_experiment = 1; % actual dose of drug 2

%% Delcare variables for synergy of output and input scores
TGI_drug1 = zeros(length(d1),length(d2)); 
TGI_drug2 = zeros(length(d1),length(d2)); 
TGI_combo = zeros(length(d1),length(d2)); 
CI_HSA = zeros(length(d1),length(d2));
CI_Bliss = zeros(length(d1),length(d2));
closest_D1 = zeros(length(d1),length(d2)); 
closest_D2 = zeros(length(d1),length(d2));
CI_LSD = zeros(length(d1),length(d2));
CI_LSD_relEC50 = zeros(length(d1),length(d2));
drug1 = zeros(length(d1),length(d2)); 
drug2 = zeros(length(d1),length(d2)); 
cum_dose = zeros(length(d1),length(d2)); 
cum_dose_rel_EC50 = zeros(length(d1),length(d2)); 
min_D = zeros(length(d1),length(d2));
min_D_relEC50 = zeros(length(d1),length(d2));
CI_Loewes = zeros(length(d1),length(d2));
onPareto = zeros(length(d1),length(d2));

%% Everything is relative to control, so set doses to 0 to find control solution
p.dose1 = 0; % mg/kg; pembro
p.dose2 = 0; % bevacizumab
sol_control  = combinations_TGI(p,in_cond,tmax,tStart1,tStart2);

%% Compute EC50 for each drug
d1_EC50 = linspace(min_d1,max_d1,num_pts*num_pts); % search over EC5-EC60 range, but a high resolution
[EC50_d1, EC50_d1_TGI] = closest_drug(d1_EC50,1,0.5,p,tmax,sol_control,in_cond,tStart1,tStart2);
fprintf('EC50 for drug 1 = %f as it has TGI = %f\n',EC50_d1, EC50_d1_TGI);
d2_EC50 = linspace(min_d2,max_d2,num_pts*num_pts); % search over EC5-EC60 range, but a high resolution
[EC50_d2, EC50_d2_TGI] = closest_drug(d2_EC50,2,0.5,p,tmax,sol_control,in_cond,tStart1,tStart2);
fprintf('EC50 for drug 2 = %f as it has TGI = %f\n',EC50_d2, EC50_d2_TGI);

%% Loop over drug doses: this will run for any model, provided the model 
%% returns tumor size at the terminal time point. Here, that output comes 
%% from combinations_TGI(); 
for i = 1:length(d1)
    p.dose1 = d1(i);  p.dose2 = 0;
    sol_drug1  = combinations_TGI(p,in_cond,tmax,tStart1,tStart2); % standalone drug1
    for j = 1:length(d2)
        p.dose1 = 0;  p.dose2 = d2(j);
        sol_drug2 = combinations_TGI(p,in_cond,tmax,tStart1,tStart2); % standalone drug2
        p.dose1 = d1(i);
        sol_combo = combinations_TGI(p,in_cond,tmax,tStart1,tStart2);
        
        %% Synergy of efficacy, where efficacy measured using TGI = percent cell kill relative to control
        TGI_drug1(i,j) = (sol_control(end)-sol_drug1(end))/sol_control(end);
%         fprintf('\tTGI(drug1) relative to control = %f\n',TGI_drug1(i,j));
        TGI_drug2(i,j) = (sol_control(end)-sol_drug2(end))/sol_control(end);
%         fprintf('\tTGI(drug2) relative to control = %f\n',TGI_drug2(i,j));
        TGI_combo(i,j) = (sol_control(end)-sol_combo(end))/sol_control(end);
%         fprintf('\tTGI(combo) relative to control = %f\n',TGI_combo(i,j));
        fprintf('Dose(drug1) = %f, Dose(drug2) = %f has TGI(combo) = %f\n',...
            p.dose1,p.dose2,TGI_combo(i,j)); 

        % HSA
        CI_HSA(i,j) = max(TGI_drug1(i,j),TGI_drug2(i,j))/TGI_combo(i,j);
        % Bliss
        bliss = TGI_drug1(i,j) + TGI_drug2(i,j) - (TGI_drug1(i,j)*TGI_drug2(i,j));
        CI_Bliss(i,j) = bliss/TGI_combo(i,j);
               
        %% Synergy of potency
        drug1(i,j) = d1(i); 
        drug2(i,j) = d2(j); 
        cum_dose(i,j) = d1(i)+d2(j); 
        cum_dose_rel_EC50(i,j) = d1(i)/EC50_d1+d2(j)/EC50_d2; 
        
        %% Find dose of drug 1 that gives TGI_combo
        [closest_D1(i,j), TGI_closest_D1] = closest_drug(d1_search,1,TGI_combo(i,j),...
            p,tmax,sol_control,in_cond,tStart1,tStart2);        
        fprintf('\tTarget TGI = %f: Dose %f of drug 1 has TGI %f\n',...
            TGI_combo(i,j),closest_D1(i,j),TGI_closest_D1);
        
        %% Find dose of drug 2 that gives TGI_combo
        [closest_D2(i,j), TGI_closest_D2] = closest_drug(d2_search,2,TGI_combo(i,j),...
            p,tmax,sol_control,in_cond,tStart1,tStart2); 
        fprintf('\tTarget TGI = %f: Dose %f of drug 2 has TGI %f\n',...
            TGI_combo(i,j),closest_D2(i,j),TGI_closest_D2); 
        
        % LSD
        min_D(i,j) = min(closest_D1(i,j),closest_D2(i,j));
        CI_LSD(i,j) = cum_dose(i,j)/min_D(i,j);
        % Loewes
        CI_Loewes(i,j) = (d1(i)/closest_D1(i,j)) + (d2(j)/closest_D2(i,j));  
        % Relative LSD
        min_D_relEC50(i,j)= min(closest_D1(i,j)/EC50_d1, closest_D2(i,j)/EC50_d2);
        CI_LSD_relEC50(i,j) = cum_dose_rel_EC50(i,j)/min_D_relEC50(i,j);
    end
end

%% TGI for drug 1
figure; 
imagesc(d1,d2,TGI_drug1'); 
colorbar
caxis([0,1]);
xlabel('Dose pembrolizumab (mg/kg)','FontSize',16); % d1
ylabel('Dose bevacizumab (mg/kg)','FontSize',16); % d2
grid off;
title('TGI(drug1) Relative to Control','FontSize',16)

%% TGI for drug 2
figure; 
imagesc(d1,d2,TGI_drug2'); 
colorbar
caxis([0,1]);
xlabel('Dose pembrolizumab (mg/kg)','FontSize',16); % d1
ylabel('Dose bevacizumab (mg/kg)','FontSize',16); % d2
grid off;
title('TGI(drug2) Relative to Control','FontSize',16)

%% Plot TGI for combination therapy
figure; hold on;
imagesc(d1,d2,TGI_combo');  hold off
xlim([0.9*min(d1),max(d1)+0.1*min(d1)]); 
ylim([0.9*min(d2),max(d2)+0.1*min(d2)]);   
colorbar
caxis([0,1]);
xlabel('Dose pembrolizumab (mg/kg)','FontSize',16); % d1
ylabel('Dose bevacizumab (mg/kg)','FontSize',16); % d2
grid off;
title('Combination TGI Relative to Control','FontSize',16)
fname_fig = [path '/TGI']; 
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);

%% Plot CI(HSA): synergy of efficacy score
% Use consistent colorbar across synergy of efficacy heatmaps
HSA_min = min(CI_HSA,[],'all');
HSA_max = max(CI_HSA,[],'all');
Bliss_min = min(CI_Bliss,[],'all'); 
Bliss_max = max(CI_Bliss,[],'all');
CI_eff_min = min(HSA_min,Bliss_min);
CI_eff_max = max(HSA_max,Bliss_max);
figure; hold on;
imagesc(d1,d2,CI_HSA'); hold off
xlim([0.9*min(d1),max(d1)+0.1*min(d1)]); 
ylim([0.9*min(d2),max(d2)+0.1*min(d2)]);   
colorbar
caxis([CI_eff_min,CI_eff_max]);
xlabel('Dose pembrolizumab (mg/kg)','FontSize',16); % d1
ylabel('Dose bevacizumab (mg/kg)','FontSize',16); % d2
grid off;
title('Synergy of Efficacy: CI(HSA)','FontSize',16)
fname_fig = [path '/syn_efficacy_CI_HSA'];
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);

%% Plot CI(Bliss): synergy of efficacy score
figure; hold on;
imagesc(d1,d2,CI_Bliss'); hold off
xlim([0.9*min(d1),max(d1)+0.1*min(d1)]); 
ylim([0.9*min(d2),max(d2)+0.1*min(d2)]);  
colorbar
caxis([CI_eff_min,CI_eff_max]);
xlabel('Dose pembrolizumab (mg/kg)','FontSize',16); % d1
ylabel('Dose bevacizumab (mg/kg)','FontSize',16); % d2
grid off;
title('Synergy of Efficacy: CI(Bliss)','FontSize',16)
fname_fig = [path '/syn_efficacy_CI_Bliss'];
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);

%% Plot cumulative drug dose relative to EC50
figure; hold on;
imagesc(d1,d2,cum_dose_rel_EC50'); hold off
xlim([0.9*min(d1),max(d1)+0.1*min(d1)]); 
ylim([0.9*min(d2),max(d2)+0.1*min(d2)]);   
colorbar
xlabel('Dose pembrolizumab (mg/kg)','FontSize',16); % d1
ylabel('Dose bevacizumab (mg/kg)','FontSize',16); % d2
grid off;
title('Cumulative Drug Dose Relative to EC50','FontSize',16)
fname_fig = [path '/cumulative_dose_relative'];
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);


%% Plot CI(LSD rel to EC50): synergy of potency score
% Use consistent colorbar across synergy of potency heatmaps
LDSrel_min = min(CI_LSD_relEC50,[],'all');
LDSrel_max = max(CI_LSD_relEC50,[],'all');
Loewe_min = min(CI_Loewes,[],'all'); 
Loewe_max = max(CI_Loewes,[],'all');
CI_pot_min = min(LDSrel_min,Loewe_min);
CI_pot_max = max(LDSrel_max,Loewe_max);
figure; hold on;
imagesc(d1,d2,CI_LSD_relEC50'); hold off
xlim([0.9*min(d1),max(d1)+0.1*min(d1)]); 
ylim([0.9*min(d2),max(d2)+0.1*min(d2)]);   
colorbar
caxis([CI_pot_min,CI_pot_max]);
xlabel('Dose pembrolizumab (mg/kg)','FontSize',16); % d1
ylabel('Dose bevacizumab (mg/kg)','FontSize',16); % d2
grid off;
title('Synergy of Potency: CI(LSD)','FontSize',16)
fname_fig = [path '/syn_potency_CI_LSD'];
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);

%% Plot CI(Loewe's): synergy of potency score
figure; hold on;
imagesc(d1,d2,CI_Loewes'); hold off
xlim([0.9*min(d1),max(d1)+0.1*min(d1)]); 
ylim([0.9*min(d2),max(d2)+0.1*min(d2)]);   
colorbar
caxis([CI_pot_min,CI_pot_max]);
xlabel('Dose pembrolizumab (mg/kg)','FontSize',16); % d1
ylabel('Dose bevacizumab (mg/kg)','FontSize',16); % d2
grid off;
title('Synergy of Potency: CI(Loewe)','FontSize',16)
fname_fig = [path '/syn_potency_CI_Loewes'];
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);

%% Isobolograms
[X,Y] = meshgrid(d1,d2);
TGI_target_vec = 0.1:0.1:0.9; 
figure;
[M,c]=contourf(X,Y,TGI_combo',TGI_target_vec,'ShowText','on');
c.LineWidth = 3;
clabel(M,c,'FontSize',15);
colorbar
xlabel('Dose pembrolizumab (mg/kg)','FontSize',16); % d1
ylabel('Dose bevacizumab (mg/kg)','FontSize',16); % d2
title('Isobologram','FontSize',16); 
fname_fig = [path '/isobologram'];
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);
 
%% Pareto front in LSD-HSA multi-synergy criterion space
CI_LSD_vec = reshape(CI_LSD_relEC50',1,[]); % Convert matrix to vector data
CI_HSA_vec = reshape(CI_HSA',1,[]);
fprintf('Pareto Front in LSDrel-HSA Multi-Synergy Criterion Space:\n'); 
[pareto_HSALSD_LSD_plot,pareto_HSALSD_HSA_plot,onPareto_LSD_HSA,onPareto] = ...
    find_Pareto(CI_LSD_vec,CI_HSA_vec,CI_LSD_relEC50,CI_HSA,drug1,drug2,onPareto);
figure; % plot objective function value at each drug dose tested
plot(CI_LSD_vec,CI_HSA_vec,'ob','MarkerFaceColor','b');  
hold on; 
plot(pareto_HSALSD_LSD_plot,pareto_HSALSD_HSA_plot,'-or','MarkerFaceColor','r'); hold off;
xlabel('Synergy of potency: CI(LSD)','FontSize',16); % objective function 2
ylabel('Synergy of efficacy: CI(HSA)','FontSize',16); % objective function 1
title('Pareto Front (Red)','FontSize',16);
fname_fig = [path '/pareto_LSD_HSA'];
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);

%% Doses on Pareto front in LSD-HSA multi-synergy criterion space
figure; hold on;
imagesc(d1,d2,TGI_combo'); 
colorbar
caxis([0,1]);
xlabel('Dose pembrolizumab (mg/kg)','FontSize',16); % d1
ylabel('Dose bevacizumab (mg/kg)','FontSize',16); % d2
grid off;
title('x = on LSD-HSA Pareto, Color=TGI','FontSize',16);
for i = 1:length(d1)
    for j = 1:length(d2)
        if(onPareto_LSD_HSA(i,j)>0)
%             fprintf('\tDose (d1,d2)=(%f,%f) appears on Pareto front\n',...
%                 drug1(i,j),drug2(i,j));
            plot(drug1(i,j),drug2(i,j),'xk','linewidth',2);
        end
    end
end
plot(d1_experiment,d2_experiment,'*r','linewidth',0.5);
hold off; 
xlim([0.9*min(d1),max(d1)+0.1*min(d1)]); 
ylim([0.9*min(d2),max(d2)+0.1*min(d2)]); 
fname_fig = [path '/doses_onPareto_LSD_HSA'];
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);


%% Pareto front in LSD-Bliss multi-synergy criterion space
CI_Bliss_vec = reshape(CI_Bliss',1,[]); % Convert matrix to vector data
fprintf('Pareto Front in LSD-Bliss Multi-Synergy Criterion Space:\n'); 
[pareto_BlissLSD_LSD_plot,pareto_BlissLSD_Bliss_plot,onPareto_LSD_Bliss,onPareto] = ...
    find_Pareto(CI_LSD_vec,CI_Bliss_vec,CI_LSD_relEC50,CI_Bliss,drug1,drug2,onPareto);
figure; % plot objective function value at each drug dose tested
plot(CI_LSD_vec,CI_Bliss_vec,'ob','MarkerFaceColor','b');  
hold on; 
plot(pareto_BlissLSD_LSD_plot,pareto_BlissLSD_Bliss_plot,'-or','MarkerFaceColor','r'); hold off;
xlabel('Synergy of potency: CI(LSD)','FontSize',16); % objective function 2
ylabel('Synergy of efficacy: CI(Bliss)','FontSize',16); % objective function 1
title('Pareto Front (Red)','FontSize',16);
fname_fig = [path '/pareto_LSD_Bliss'];
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);

%% Doses on Pareto front in LSD-Bliss multi-synergy criterion space
figure; hold on;
imagesc(d1,d2,TGI_combo'); 
colorbar
caxis([0,1]);
xlabel('Dose pembrolizumab (mg/kg)','FontSize',16); % d1
ylabel('Dose bevacizumab (mg/kg)','FontSize',16); % d2
grid off;
title('x = on LSD-Bliss Pareto, Color=TGI','FontSize',16);
for i = 1:length(d1)
    for j = 1:length(d2)
        if(onPareto_LSD_Bliss(i,j)>0)
%             fprintf('\tDose (d1,d2)=(%f,%f) appears on Pareto front\n',...
%                 drug1(i,j),drug2(i,j));
            plot(drug1(i,j),drug2(i,j),'xk','linewidth',2);
        end
    end
end
plot(d1_experiment,d2_experiment,'*r','linewidth',0.5);
hold off; 
xlim([0.9*min(d1),max(d1)+0.1*min(d1)]); 
ylim([0.9*min(d2),max(d2)+0.1*min(d2)]); 
fname_fig = [path '/doses_onPareto_LSD_Bliss'];
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);

%% Pareto front in Loewes-HSA multi-synergy criterion space
CI_Loewes_vec = reshape(CI_Loewes',1,[]); % Convert matrix to vector data
fprintf('Pareto Front in Loewes-HSA Multi-Synergy Criterion Space:\n'); 
[pareto_LoewesHSA_Loewes_plot,pareto_LoewesHSA_HSA_plot,onPareto_Loewes_HSA,onPareto] = ...
    find_Pareto(CI_Loewes_vec,CI_HSA_vec,CI_Loewes,CI_HSA,drug1,drug2,onPareto);
figure; % plot objective function value at each drug dose tested
plot(CI_Loewes_vec,CI_HSA_vec,'ob','MarkerFaceColor','b');  
hold on; 
plot(pareto_LoewesHSA_Loewes_plot,pareto_LoewesHSA_HSA_plot,'-or','MarkerFaceColor','r'); hold off;
xlabel('Synergy of potency: CI(Loewe)','FontSize',16); % objective function 2
ylabel('Synergy of efficacy: CI(HSA)','FontSize',16); % objective function 1
title('Pareto Front (Red)','FontSize',16);
fname_fig = [path '/pareto_Loewes_HSA'];
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);

%% Doses on Pareto front in Loewes-HSA multi-synergy criterion space
figure; hold on;
imagesc(d1,d2,TGI_combo'); 
colorbar
caxis([0,1]);
xlabel('Dose pembrolizumab (mg/kg)','FontSize',16); % d1
ylabel('Dose bevacizumab (mg/kg)','FontSize',16); % d2
grid off;
title('x = on Loewe-HSA Pareto, Color=TGI','FontSize',16);
for i = 1:length(d1)
    for j = 1:length(d2)
        if(onPareto_Loewes_HSA(i,j)>0)
            fprintf('\tDose (d1,d2)=(%f,%f) appears on Pareto front\n',...
                drug1(i,j),drug2(i,j));
            plot(drug1(i,j),drug2(i,j),'xk','linewidth',2);
        end
    end
end
plot(d1_experiment,d2_experiment,'*r','linewidth',0.5);
hold off; 
xlim([0.9*min(d1),max(d1)+0.1*min(d1)]); 
ylim([0.9*min(d2),max(d2)+0.1*min(d2)]); 
fname_fig = [path '/doses_onPareto_Loewes_HSA'];
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);

%% Pareto front in Loewes-Bliss multi-synergy criterion space
fprintf('Pareto Front in Loewes-Bliss Multi-Synergy Criterion Space:\n'); 
[pareto_LoewesBliss_Loewes_plot,pareto_LoewesBliss_Bliss_plot,onPareto_Loewes_Bliss,onPareto] = ...
    find_Pareto(CI_Loewes_vec,CI_Bliss_vec,CI_Loewes,CI_Bliss,drug1,drug2,onPareto);
figure; % plot objective function value at each drug dose tested
plot(CI_Loewes_vec,CI_Bliss_vec,'ob','MarkerFaceColor','b');  
hold on; 
plot(pareto_LoewesBliss_Loewes_plot,pareto_LoewesBliss_Bliss_plot,'-or','MarkerFaceColor','r'); hold off;
xlabel('Synergy of potency: CI(Loewe)','FontSize',16); % objective function 2
ylabel('Synergy of efficacy: CI(Bliss)','FontSize',16); % objective function 1
title('Pareto Front (Red)','FontSize',16);
fname_fig = [path '/pareto_Loewes_Bliss'];
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);

%% Doses on Pareto front in Loewes-Bliss multi-synergy criterion space
figure; hold on;
imagesc(d1,d2,TGI_combo'); 
colorbar
caxis([0,1]);
xlabel('Dose pembrolizumab (mg/kg)','FontSize',16); % d1
ylabel('Dose bevacizumab (mg/kg)','FontSize',16); % d2
grid off;
title('x = on Loewe-Bliss Pareto, Color=TGI','FontSize',16);
for i = 1:length(d1)
    for j = 1:length(d2)
        if(onPareto_Loewes_Bliss(i,j)>0)
            fprintf('\tDose (d1,d2)=(%f,%f) appears on Pareto front\n',...
                drug1(i,j),drug2(i,j));
            plot(drug1(i,j),drug2(i,j),'xk','linewidth',2);
        end
    end
end
plot(d1_experiment,d2_experiment,'*r','linewidth',0.5);
hold off; 
xlim([0.9*min(d1),max(d1)+0.1*min(d1)]); 
ylim([0.9*min(d2),max(d2)+0.1*min(d2)]); 
fname_fig = [path '/doses_onPareto_Loewes_Bliss'];
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);

%% Doses that appear on any Pareto front
count_onPareto = zeros(length(d1),length(d2)); 
for i = 1:length(d1)
    for j = 1:length(d2)
        if(onPareto_LSD_HSA(i,j)>0)
            count_onPareto(i,j) = count_onPareto(i,j) + 1; 
        end
        if(onPareto_LSD_Bliss(i,j)>0)
            count_onPareto(i,j) = count_onPareto(i,j) + 1;
        end
        if(onPareto_Loewes_HSA(i,j)>0)
            count_onPareto(i,j) = count_onPareto(i,j) + 1; 
        end
        if(onPareto_Loewes_Bliss(i,j)>0)
            count_onPareto(i,j) = count_onPareto(i,j) + 1; 
        end
    end
end

font = 'Arial';
step_d1 = 0.5*(max(d1)-min(d1))/num_pts; 
step_d2 = 0.5*(max(d2)-min(d2))/num_pts; 
figure; hold on;
imagesc(d1,d2,TGI_combo'); 
xlim([min(d1)-step_d1,max(d1)+step_d1]); 
ylim([min(d2)-step_d2,max(d2)+step_d2]); 
colorbar
caxis([0,1]);
xlabel('Dose pembrolizumab (mg/kg)','FontSize',16); % d1
ylabel('Dose bevacizumab (mg/kg)','FontSize',16); % d2
grid off;
title('Number of Times Dose Appears on Pareto Front','FontSize',16);
subtitle('Color Gives TGI of Combination','FontSize',14);
for i = 1:length(d1)
    for j = 1:length(d2)
       m = num2str(count_onPareto(i,j));
       if(count_onPareto(i,j) > 0)
           text(drug1(i,j),drug2(i,j),m,'fontname',font,'color','black','fontweight','bold');
       end
    end
end
plot(d1_experiment,d2_experiment,'*r','linewidth',0.5);
hold off; 
fname_fig = [path '/doses_onPareto_any'];
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);


%% Display all points that appear on multiple Pareto fronts
fprintf('\nDoses that appear on multiple Pareto fronts:\n');
for i = 1:length(d1)
    for j = 1:length(d2)
        if(onPareto(i,j)>1)
            fprintf('\tDose (d1,d2)=(%f,%f) appears on %d Pareto fronts:',...
                drug1(i,j),drug2(i,j),onPareto(i,j));
            if(onPareto_LSD_HSA(i,j)==1)
                fprintf(' LSD-HSA,');
            end
            if(onPareto_LSD_Bliss(i,j)==1)
                fprintf(' LSD-Bliss,');
            end
            if(onPareto_Loewes_HSA(i,j)==1)
                fprintf(' Loewes-HSA,');
            end
            if(onPareto_Loewes_Bliss(i,j)==1)
                fprintf(' Loewes-Bliss,');
            end
            fprintf('\n');     
        end
    end
end

dose_exp = [d1_experiment d2_experiment];
distPareto = 10^10; d1_idx = -1; d2_idx = -1;  
dose_exp = [d1_experiment d2_experiment];
distPareto = 10^10; d1_idx = -1; d2_idx = -1;  
for i = 1:length(d1)
    for j = 1:length(d2)
        if(count_onPareto(i,j) > 0)
            dose_diff = dose_exp-[d1(i) d2(j)];
            dtemp = norm(dose_diff);
            if dtemp<distPareto
               distPareto = dtemp;
               d1_idx = d1(i);
               d2_idx = d2(j);
            end
        end
    end
end
fprintf('The Pareto optimal dose closest to experimental dose gives d1 = %f, d2 = %f\n',...
    d1_idx,d2_idx); 

toc
diary off
fsave = [path '/output.mat']; 
save(fsave,'TGI_drug1','TGI_drug2','TGI_combo','CI_HSA','CI_Bliss','closest_D1', ...
    'closest_D2','CI_LSD','CI_Loewes','drug1','drug2','cum_dose','min_D', ...
    'pareto_HSALSD_LSD_plot','pareto_HSALSD_HSA_plot','pareto_BlissLSD_LSD_plot',...
    'pareto_BlissLSD_Bliss_plot','pareto_LoewesHSA_Loewes_plot', ...
    'pareto_LoewesHSA_HSA_plot','pareto_LoewesBliss_Loewes_plot',...
    'pareto_LoewesBliss_Bliss_plot','onPareto','onPareto_LSD_HSA',...
    'onPareto_LSD_Bliss','onPareto_Loewes_HSA','onPareto_Loewes_Bliss', ...
    'CI_LSD_relEC50','min_D_relEC50','d1_idx','d2_idx','cum_dose_rel_EC50'); 
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions needed to run script                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define model parameters, initial conditions, and dosing protocol
%% **User must set these for their particular model **

function [p,in_cond,tmax,tStart1,tStart2,d1,d2,d1_search,d2_search] = ...
    set_parameters_ICs_protocol(d1_spacing,d2_spacing,num_pts)
    tmax = 35; 
    % drug 1 - pembro, 10mpk
    p.dose1 = 10; % mg/kg; pembro
    p.dosenumD1 = 5;   % number of doses
    p.freq1=1; % day
    p.intvlD1 = d1_spacing*p.freq1; % dosing interval, days
    tStart1 = 3; % drug 1, time of start
    % drug 2 - bevacizumab, 1mpk
    p.dose2 = 1; % bevacizumab
    p.dosenumD2 = 6;   % number of doses
    p.freq2=1; % day
    p.intvlD2 = d2_spacing*p.freq2; % dosing interval, days
    tStart2 = 0.00001; % drug 2, time of start
    p.iv = 0; % switch for IV (p.iv==1) vs SC (p.iv ==0)
    % switch for cell line, 1 = H1299 (better fits), 2 = A549 (worse fits)
    p.cell_line = 2;
%     if(p.cell_line == 1)
%         fprintf('Cell line %d = H1299\n',p.cell_line);
%     else
%         fprintf('Cell line %d = A549\n',p.cell_line);
%     end

    %% Initial conditions
    D1_iv0 = 0; D1_p0  = 0; D1_t0  = 0; % drug 1
    D2_iv0 = 0; D2_p0  = 0; 
    x20    = 0; x30    = 0; x40    = 0; %"transitional tumor volume things"
    if p.cell_line == 1
        x10 = 36.86513265; %initial tumor volume
    elseif p.cell_line == 2
        x10 = 43.62606232;
    end
    % these are currently irrelevant, free variables for potential future modifications
    D2_t0  = 0; D2T20  = 0; %drug 2
    D1T10  = 0; %drug 1
    T10    = 0; %concentration of target 1
    T20    = 0; %concentration of target 2
    in_cond=[D1_iv0 D1_p0 D1_t0 D1T10 T10 D2_iv0 D2_p0 D2_t0 D2T20 T20 x10 x20 x30 x40];

    %% Parameters
    %% PK parameters
    % drug 1 pembro, 2 compartment
    p.V1_1  = 70; % ml/kg 
    p.V2_1  = 33;% ml/kg
    if p.iv == 1
        p.k01_1 = 0;
    elseif p.iv == 0
        p.k01_1 = 0.11; % SC
    end
    p.Cl1_1 = 20;% ml/kg/day
    p.Cl2_1 = 22;% ml/kg/day
    p.k10_1 = p.Cl1_1/p.V1_1;
    p.k12_1 = p.Cl2_1/p.V1_1;
    p.k21_1 = p.Cl2_1/p.V2_1;

    % drug 2, bevacizumab, 1 compartment
    p.V1_2  = 119; % ml/kg
    p.Cl1_2 = 13.6 ; % ml/kd/day
    p.k10_2 = p.Cl1_2/p.V1_2;
    if p.iv == 1
        p.k01_2 = 0;
    elseif p.iv == 0
        p.k01_2 = 1.3; % SC
    end

    %% PD parameters
    if p.cell_line == 1
        p.lamba1 = .11;
    elseif p.cell_line == 2
        p.lamba1 = .11;
    end
    p.K = 100000;

    %% Combination efficacy parameters
    if p.cell_line == 1
        p.k2_d1 = 0.0008; % effect of pembro 10mpk
        p.k2_d2 = 0.0011; % effect of beva 1 mpk
        p.k3 = 0.0001; %combination factor
    elseif p.cell_line == 2 %A549
        p.k2_d1 = 0.0001; % effect of pembro 10mpk
        p.k2_d2 = 0.0008; % effect of beva 1 mpk
        p.k3 = 0.00002;
    end

    % Parameters for potential modification of relative contribution of drugs
    p.alpha = 1;
    p.beta  = 1;
    p.k1 = .000575; %transitional states
    
    %% Finally, set d1 and d2 range to search
    if(p.cell_line == 1) % H1299
        min_d1 = 0.408; max_d1 = 10.88; % Experiments use 10 mg/kg of pembro (up to EC70)
        d1 = linspace(min_d1,max_d1,num_pts); % EC5 to EC6x for drug 1 monotherapy
        min_d2 = 0.1372; max_d2 = 3.457; % Experiments use 1 mpk bevacizumab (up to EC70) 
        d2 = linspace(min_d2,max_d2,num_pts); % EC5 to EC6x for drug 2 monotherapy
        %% From EC70 range: E(max_d1,max_d2)=0.9495, which is achieved if go up to monotherapy D1=120, D2=15
        % Dose range for drug 1 that ensures we can get monotherapy efficacy = E(max_d1,max_d2)
        d1_search = linspace(min(d1),120,30*num_pts); 
        % Dose range for drug 2 that ensures we can get monotherapy efficacy = E(max_d1,max_d2)  
        d2_search = linspace(min(d2),15,30*num_pts); 
    else % A549
        min_d1 = 3.275; max_d1 = 87.2; % Experiments use 10 mg/kg of pembro
        d1 = linspace(min_d1,max_d1,num_pts); % EC5 to EC60 for drug 1 monotherapy
        min_d2 = 0.1892; max_d2 = 4.762; % Experiments use 1 mpk bevacizumab 
        d2 = linspace(min_d2,max_d2,num_pts); % EC5 to EC60 for drug 2 monotherapy
        %% E(max_d1,max_d2)=0.9559, which is achieved if go up to monotherapy D1=1650, D2=20
        % Dose range for drug 1 that ensures we can get monotherapy efficacy = E(max_d1,max_d2)
        d1_search = linspace(min(d1),1650,30*num_pts); 
        % Dose range for drug 2 that ensures we can get monotherapy efficacy = E(max_d1,max_d2)  
        d2_search = linspace(min(d2),20,30*num_pts); 
    end
    
end

%% Find dose of specified monotherapy closest to target TGI
function [closest_D, closest_D_TGI] = closest_drug(d,d_num,TGI_target,p,...
    tmax,sol_control,in_cond,tStart1,tStart2)
    dist_closest_D = 10^10; closest_D = -1; closest_D_TGI = -1; 
    for k = 1:length(d)
        if d_num==1            
            p.dose1 = d(k);  p.dose2 = 0;
            sol_d = combinations_TGI(p,in_cond,tmax,tStart1,tStart2); % standalone drug1
            
        elseif d_num==2
            p.dose1 = 0; p.dose2 = d(k);  
            sol_d = combinations_TGI(p,in_cond,tmax,tStart1,tStart2); % standalone drug1
        else
            fprintf('Can only use d_num = 1 or 2, but set d_num = %d\n',d_num);
            STOP
        end
        TGI_d = (sol_control(end)-sol_d(end))/sol_control(end);
        dist_d = abs(TGI_d-TGI_target); 
        if(dist_d<dist_closest_D)
            closest_D = d(k); 
            closest_D_TGI = TGI_d; 
            dist_closest_D = dist_d; 
        end
    end
end

%% Find points on Pareto front
function [par_potency_plot,par_efficacy_plot,pareto,onPareto] = ...
    find_Pareto(CI_potency_vec,CI_efficacy_vec,CI_pot,CI_eff,drug1,drug2,onPareto)

    %% Find entire boundary
    w = boundary(CI_potency_vec',CI_efficacy_vec',1);
    
    %% Not all boundary points on Pareto front: need to minimize both objectives 
    pareto_pot = []; pareto_eff = [];
    count_par = 0;
    xbnd = CI_potency_vec(w); % reduces set to boundary points
    ybnd = CI_efficacy_vec(w);
    for q=1:length(w)
        xdist = xbnd-xbnd(q); 
        xf = find(xdist<=0); % <=0 so only look at points with same or smaller objective
        if isempty(xf) == 1 % empty so largest x point
            % fprintf('%d: (x,y)=(%f,%f): no points have smaller x\n',q,xbnd(q),ybnd(q)); 
            count_par = count_par+1; 
            pareto_pot(count_par) =  xbnd(q);
            pareto_eff(count_par) =  ybnd(q);
            % fprintf('\tAdded to Pareto front\n');
        else % not empty so now we have to ask - do any of these points have smaller y? 
            % fprintf('%d: (x,y)=(%f,%f): has points with larger x\n',q,xbnd(q),ybnd(q)); 
            ybnd2 = ybnd(xf); % reduce set further to only look at points with <= first objective
            ydist = ybnd2-ybnd(q);
            if isempty(find(ydist<0, 1)) == 1 % no points with smamler first objective have smaller second objective
                % fprintf('\tAND no points with smaller y\n');
                count_par = count_par+1; 
                pareto_pot(count_par) =  xbnd(q);
                pareto_eff(count_par) =  ybnd(q);
                % fprintf('\tAdded to Pareto front\n');
            else
                % fprintf('\tBUT points with smaller synergy of efficacy - NOT on Pareto front\n');
            end 
        end
    end
    % Sort points in order of increasing x (first objective) for connecting plot
    [par_potency_plot, index] = sort(pareto_pot); 
    par_efficacy_plot = pareto_eff(index); 

    %% Save and list all points on Pareto front
    pareto = zeros(length(drug1),length(drug2));
    for k = 1:length(pareto_pot)
        fprintf('\tPareto point #%d:\n',k); 
        [I_pot_row,I_pot_col] = find(CI_pot==pareto_pot(k));
        [I_eff_row,I_eff_col] = find(CI_eff==pareto_eff(k));
        Iintersect_row = intersect(I_pot_row,I_eff_row);
        Iintersect_col = intersect(I_pot_col,I_eff_col);
        if (numel(Iintersect_row)==1)&&(numel(Iintersect_col)==1)
            fprintf('\t\t"Dose": Drug 1 = %f, Drug 2 = %f\n',...
                drug1(Iintersect_row,Iintersect_col),drug2(Iintersect_row,Iintersect_col));
            fprintf('\t\tScores: CI(potency) = %f, CI(efficacy) = %f\n',...
                CI_pot(Iintersect_row,Iintersect_col),CI_eff(Iintersect_row,Iintersect_col));
            onPareto(Iintersect_row,Iintersect_col) = onPareto(Iintersect_row,Iintersect_col) + 1;
            pareto(Iintersect_row,Iintersect_col) = 1; 
        else 
            fprintf('Did NOT find a unique protocol corresponding to CI(potency) = %f and CI(efficacy) = %f\n',...
                pareto_pot(k),pareto_eff(k));
            STOP
        end
    end
end

%% PK/PD model of pembro+beva. Just set to return tumor size at end of tspan
%% **User must set rewrite this for their particular model **
function tumor = combinations_TGI(p,y0,tspan,tStart1,tStart2)
    function dydt = ode(t,y0)
        D1_iv = y0(1);
        D1_p  = y0(2);
        D1_t  = y0(3);
        D1T1  = y0(4);
        T1    = y0(5);
        
        D2_iv = y0(6);
        D2_p  = y0(7);
        D2_t  = y0(8);
        D2T2  = y0(9);
        T2    = y0(10);
        
        x1    = y0(11);
        x2    = y0(12);
        x3    = y0(13);
        x4    = y0(14);
        wt=x1+x2+x3+x4;
        
        
        %% drug 1
        % pembro, 2 compartment       
        dD1_iv=-p.k01_1*D1_iv;
        dD1_p=p.k01_1*D1_iv-p.k10_1*D1_p...
            -p.k12_1*D1_p+p.k21_1*D1_t*p.V2_1/p.V1_1;
        dD1_t=p.k12_1*D1_p*(p.V1_1/p.V2_1) - p.k21_1*D1_t;
        dT1D1 = 1; %free, to have option to add interaction with target
        dT1 = 1; % %free, to have option to add target turnover
        
        %% drug 2, bevacizumab, 1 compartment
        dD2_iv=-p.k01_2*D2_iv; %absorption
        dD2_p=p.k01_2*D2_iv-p.k10_2*D2_p; %distribution and clearance
        dD2_t=1; %free, to have option to do 2 compartment        
        dT2D2 = 1; %free, to have option to add target turnover        
        dT2 = 1; %        
        kill_term = p.k2_d1*D1_p + p.k2_d2*D2_p...
            + p.k3*(D1_p^p.alpha)*(D2_p^p.beta);
        %% tumor
        %%%%%%%%% the important equation %%%%%%%
        %    if p.case == 1; if logistic, get ...%
        dx1 = p.lamba1*x1*(1-wt/p.K) - x1*kill_term;...%p.lamda0*x1/((1.0+(p.lamda0/p.lamda1*wt)^p.psi)^(1.0/p.psi))... %typical simeoni, what we usually use
        %             - p.k2_d1*D1_p*x1 - p.k2_d2*D2_p*x1...
        %             - p.k3*(D1_p^p.alpha)*(D2_p^p.beta)*x1; %just regular "additive"
        %   end
        %%%%%%%%% the important equation %%%%%%%
        dx2 = x1*kill_term - p.k1*x2;
        dx3 = p.k1*(x2 -x3);
        dx4 = p.k1*(x3 -x4);
             
        %%
        dydt=[dD1_iv;dD1_p; dD1_t; dT1D1; dT1;... % drug 1 and target 1
            dD2_iv;dD2_p; dD2_t; dT2D2; dT2;... % drug 2 and target 2
            dx1 ; dx2 ; dx3; dx4]; % tumor
    end

    % dose recalculations fron to ng/ml here because need V1 defined
    p.doseD1 = p.dose1*10^3/p.V1_1;%/p.MW1;% mpk, 
    p.doseD2 = p.dose2*10^3/p.V1_2;%/p.MW2;% mpk, 

    ivdoseD1 = p.doseD1;
    ivdoseD2 = p.doseD2;

    D1_iv0 = y0(1);
    D1_p0  = y0(2);
    D1_t0  = y0(3);
    D1T10  = y0(4);
    T10    = y0(5);

    D2_iv0 = y0(6);
    D2_p0  = y0(7);
    D2_t0  = y0(8);
    D2T20  = y0(9);
    T20    = y0(10);

    x10    = y0(11);
    x20    = y0(12);
    x30    = y0(13);
    x40    = y0(14);
    % 
    for i = 1:p.dosenumD1
        tDose1(i) = p.intvlD1*(i-1) + tStart1;
        cDose1(i) = ivdoseD1;
    end

    % If you change p.dosenumD2 to something other than
    % p.dosenumD1, then the lengths of tDose2 and cDose2
    % are different from tDose1 and cDose1

    %     tDose2(1) = tStart2;
    %     cDose2(1) = p.doseD2_load;
    % for i = 2:p.dosenumD2
    for i = 1:p.dosenumD2
        tDose2(i) = p.intvlD2*(i-1) + tStart2;
        cDose2(i) = ivdoseD2;
    end

    % Your time interval is the UNION of both tDose1 and tDose2
    % Put the 2 sets together, in order, no duplicates
    tUnion = union( tDose1, tDose2, 'sorted');
    tUnion = [0 tUnion];
    tLength = length(tUnion);
    ti = [tUnion max(tspan)];

    % Run the for the union of time; step through each interval
    %for i = 1:length(tDose1)
    tPoints = []; cPoints = [];
    for i = 1:tLength

        %cc_0 = cc_0 + cDose(i);
        % Given a TIME, does that TIME exist in tDose1 array?
        [lia1,locb1] = ismember( tUnion(i), tDose1 );
        if lia1 == 1
            if p.iv == 1 %IV
                D1_p0 = D1_p0 + cDose1(locb1); %give drug 1 (Orencia)
            elseif p.iv == 0 %SC
                D1_iv0 = D1_iv0 + cDose1(locb1); %give drug 1 (Orencia)
            end
        end

        % Does that TIME exist in tDose2 array?
        [lia2,locb2] = ismember( tUnion(i), tDose2 );
        if lia2 == 1
            if p.iv == 1 %IV
                D2_p0 = D2_p0 + cDose2(locb2); %give drug 2 
            elseif p.iv == 0
                D2_iv0 = D2_iv0 + cDose2(locb2); %give drug 2 
            end
        end


        % Build the tSpan and run the Diff Eq
        tSpan = [ti(i) ti(i+1)];
        [tOut,cOut] = ode45(@ode,tSpan,[D1_iv0 D1_p0 D1_t0 D1T10 T10 D2_iv0 D2_p0 D2_t0 D2T20 T20 x10 x20 x30 x40]);

        % Concatenate the results of the ode45 work
        tPoints = [tPoints; tOut];
        cPoints = [cPoints; cOut];


        % Set new init conditions for ode45 and repeat the loop
        D1_iv0 = cOut(end,1); D1_p0 = cOut(end,2); D1_t0 = cOut(end,3);...
            D1T10 = cOut(end,4); T10 = cOut(end,5);...
            D2_iv0 = cOut(end,6); D2_p0 = cOut(end,7); D2_t0 = cOut(end,8);...
            D2T20 = cOut(end,9); T20 = cOut(end,10);...
            x10 = cOut(end,11); x20 = cOut(end,12);...
            x30 = cOut(end,13); x40 = cOut(end,14);
    end

    % Send these results out of the function
    X = cPoints;
    tumor = X(:,11)+X(:,12)+X(:,13)+X(:,14);
    %Time = tPoints;
end


