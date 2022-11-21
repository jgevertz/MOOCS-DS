%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% Multi-objective synergy optimization by Jana Gevertz (8/30/2022)      %
% - The code first establishes a tumor that grows logistically at rate  % 
%   r. Drug 1 kills at rate d1, drug 2 at rate d2                       %
% - At a fixed value of r, we sweep over a range of d1 and d2 values    %
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
%      weighted sum.                                                    %
%   4) Synergy of potency using Loewe's additivity formula:             %
%       CI(Loewes) = (d1/D1)+(d2/D2)                                    %
%      where CI = 1 is additive, CI<1 is synergy, CI>1 is antagonism    %
%      (Note: could use generalized Loewe's CI instead: see code        %
%      Understanding_Loewes/Loewes_GeneralizedLoewes_CI.m               %
% - The code plots:                                                     %        
%   Figure 1: heatmap of combination TGI relative to control            %
%   Figure 2: heatmap of synergy of efficacy defined by HSA             %
%   Figure 3: heatmap of synergy of efficacy defined by Bliss           %
%   Figure 4: heatmap of cumulative dose                                %
%   Figure 5: heatmap of synergy of potency defined by LSD              %
%   Figure 6: heatmap of synergy of potency defined by Loewe            %
%   Figure 7: Pareto front LSD-HSA multi-synergy criterion space        %
%   Figure 8: Pareto front LSD-Bliss multi-synergy criterion space      %
%   Figure 9: Pareto front Loewe-HSA multi-synergy criterion space      %
%   Figure 10: Pareto front Loewe-Bliss multi-synergy criterion space   %
%   Figure 11: Contour plot of isobolograms (curves of constant         %
%              combination efficacy)                                    %
% - All output is written to PointsOnPareto, including a list of all    %
%   doses on each Pareto front, and a list of all doses that appear on  %
%   multiple Pareto fronts                                              %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all; tic;
diary PointsOnPareto
% DE parameters and options
global r M y0 tf
r = 0.2; % growth rate
M = 10; %Carrying capacity
y0 = 1; % IC
tf = 10; % time horizon

% Range of kill rates
num_pts = 50; %30
min_d1 = 0.00684; max_d1 = 0.2225; % EC5 to EC85 for drug 1 monotherapy
d1 = linspace(min_d1,max_d1,num_pts); % drug 1: need d1>r to be effective
min_d2 = 0.002695; max_d2 = 0.2903; % EC5 to EC85 for drug 1 monotherapy
d2 = linspace(min_d2,max_d2,num_pts); % drug 2: need d2>r to be effective

%% For synergy of output and input scores
TGI_drug1 = zeros(length(d1),length(d2)); 
TGI_drug2 = zeros(length(d1),length(d2)); 
TGI_combo = zeros(length(d1),length(d2)); 
CI_HSA = zeros(length(d1),length(d2));
CI_Bliss = zeros(length(d1),length(d2));
closest_D1 = zeros(length(d1),length(d2)); 
closest_D2 = zeros(length(d1),length(d2));
CI_LSD = zeros(length(d1),length(d2));
drug1 = zeros(length(d1),length(d2)); 
drug2 = zeros(length(d1),length(d2)); 
cum_dose = zeros(length(d1),length(d2)); 
min_D = zeros(length(d1),length(d2));
CI_Loewes = zeros(length(d1),length(d2));
cum_dose_rel_EC50 = zeros(length(d1),length(d2)); 
min_D_relEC50 = zeros(length(d1),length(d2));
CI_LSD_relEC50 = zeros(length(d1),length(d2));
global onPareto
onPareto = zeros(length(d1),length(d2));

%% Everything is relative to control
global sol_control
sol_control = ode23s(@(t,x) logistic_drug(t,x,r,M,0,0),[0,tf],y0);
d1_EC50 = linspace(0.08,0.09,10*num_pts); % manualy determined EC50 in this range
[EC50_d1, EC50_d1_TGI] = closest_drug(d1_EC50,1,0.5);
d2_EC50 = linspace(0.05,0.06,10*num_pts); % manualy determined EC50 in this range
[EC50_d2, EC50_d2_TGI] = closest_drug(d2_EC50,2,0.5);

%% Loop over drug doses (really, kill rates here)
for i = 1:length(d1)
    sol_drug1 = ode23s(@(t,x) logistic_drug(t,x,r,M,d1(i),0),[0,tf],y0); % standalone drug1
    for j = 1:length(d2)
        sol_drug2 = ode23s(@(t,x) logistic_drug(t,x,r,M,0,d2(j)),[0,tf],y0); % standalone drug2
        sol_combo = ode23s(@(t,x) logistic_drug(t,x,r,M,d1(i),d2(j)),[0,tf],y0);
        
        %% Synergy of efficacy, where efficacy measured using TGI = percent cell kill relative to control
        TGI_drug1(i,j) = (sol_control.y(end)-sol_drug1.y(end))/sol_control.y(end);
%         fprintf('\tTGI(drug1) relative to control = %f\n',TGI_drug1(i,j));
        TGI_drug2(i,j) = (sol_control.y(end)-sol_drug2.y(end))/sol_control.y(end);
%         fprintf('\tTGI(drug2) relative to control = %f\n',TGI_drug2(i,j));
        TGI_combo(i,j) = (sol_control.y(end)-sol_combo.y(end))/sol_control.y(end);
%         fprintf('\tTGI(combo) relative to control = %f\n',TGI_combo(i,j));
        fprintf('Dose(drug1) = %f, Dose(drug2) = %f has TGI(combo) = %f\n',...
            d1(i),d2(j),TGI_combo(i,j)); 
        
        % Use HSA
        CI_HSA(i,j) = max(TGI_drug1(i,j),TGI_drug2(i,j))/TGI_combo(i,j);
        % Use Bliss
        bliss = TGI_drug1(i,j) + TGI_drug2(i,j) - (TGI_drug1(i,j)*TGI_drug2(i,j));
        CI_Bliss(i,j) = bliss/TGI_combo(i,j);
               
        %% Synergy of input using LSD
        drug1(i,j) = d1(i); 
        drug2(i,j) = d2(j); 
        cum_dose(i,j) = d1(i)+d2(j); 
        cum_dose_rel_EC50(i,j) = d1(i)/EC50_d1+d2(j)/EC50_d2; 
        
        %% Find dose of drug 1 that gives TGI_combo
        d3 = linspace(min(d1),0.36,10*num_pts); % manually determined upper bound so TGI(upperbound) = TGI(best combo)
        %d3 = linspace(0,2*(max(d1)+max(d2)),10*num_pts); 
        [closest_D1(i,j), TGI_closest_D1] = closest_drug(d3,1,TGI_combo(i,j));        
        fprintf('\tTarget TGI = %f: Dose %f of drug 1 has TGI %f\n',...
            TGI_combo(i,j),closest_D1(i,j),TGI_closest_D1);
        
        %% Find dose of drug 2 that gives TGI_combo
        d4 = linspace(min(d2),1.1,10*num_pts); % manually determined upper bound so TGI(upperbound) = TGI(best combo)
        [closest_D2(i,j), TGI_closest_D2] = closest_drug(d4,2,TGI_combo(i,j)); 
        fprintf('\tTarget TGI = %f: Dose %f of drug 2 has TGI %f\n',...
            TGI_combo(i,j),closest_D2(i,j),TGI_closest_D2); 
       
        %% Synergy of potency
        min_D(i,j) = min(closest_D1(i,j),closest_D2(i,j));
        CI_LSD(i,j) = cum_dose(i,j)/min_D(i,j);
        CI_Loewes(i,j) = (d1(i)/closest_D1(i,j)) + (d2(j)/closest_D2(i,j)); 
        % Use relative LSD instead (to normalize doses)
        min_D_relEC50(i,j)= min(closest_D1(i,j)/EC50_d1, closest_D2(i,j)/EC50_d2);
        CI_LSD_relEC50(i,j) = cum_dose_rel_EC50(i,j)/min_D_relEC50(i,j);
    end
end

%% Plot TGI for combination therapy
figure; 
imagesc(d1,d2,TGI_combo'); 
colorbar
caxis([0,1]);
xlabel('Killing rate of Drug 1','FontSize',16); % d1
ylabel('Killing rate of Drug 2','FontSize',16); % d2
grid off;
title('Combination TGI Relative to Control','FontSize',16)
fname_fig = 'TGI';
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);

%% Plot CI(HSA): synergy of efficacy score
figure; 
imagesc(d1,d2,CI_HSA'); 
colorbar
xlabel('Killing rate of Drug 1','FontSize',16); % d1
ylabel('Killing rate of Drug 2','FontSize',16); % d2
grid off;
title('Synergy of Efficacy: CI(HSA)','FontSize',16)
fname_fig = 'syn_efficacy_CI_HSA';
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);

%% Plot CI(Bliss): synergy of efficacy score
figure; 
imagesc(d1,d2,CI_Bliss'); 
colorbar
xlabel('Killing rate of Drug 1','FontSize',16); % d1
ylabel('Killing rate of Drug 2','FontSize',16); % d2
grid off;
title('Synergy of Efficacy: CI(Bliss)','FontSize',16)
fname_fig = 'syn_efficacy_CI_Bliss';
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);

%% Plot cumulative drug dose
figure; 
imagesc(d1,d2,cum_dose'); 
colorbar
caxis([0,d1(end)+d2(end)]);
xlabel('Killing rate of Drug 1','FontSize',16); % d1
ylabel('Killing rate of Drug 2','FontSize',16); % d2
grid off;
title('Cumulative Drug Dose','FontSize',16)
fname_fig = 'cumulative_dose';
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);

%% Plot cumulative drug dose relative to EC50
figure; 
imagesc(d1,d2,cum_dose_rel_EC50'); 
colorbar
%caxis([0,d1(end)+d2(end)]);
xlabel('Killing rate of Drug 1','FontSize',16); % d1
ylabel('Killing rate of Drug 2','FontSize',16); % d2
grid off;
title('Cumulative Drug Dose Relative to EC50','FontSize',16)
fname_fig = 'cumulative_dose_relative';
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);

%% Plot CI(LSD): synergy of potency score
figure;
imagesc(d1,d2,CI_LSD'); 
colorbar
xlabel('Killing rate of Drug 1','FontSize',16); % d1
ylabel('Killing rate of Drug 2','FontSize',16); % d2
grid off;
title('Synergy of Potency: CI(LSD)','FontSize',16)
fname_fig = 'syn_potencyy_CI_LSD';
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);

figure;
imagesc(d1,d2,CI_LSD_relEC50'); 
colorbar
%caxis([CI_pot_min,CI_pot_max]);
xlabel('Killing rate of Drug 1','FontSize',16); % d1
ylabel('Killing rate of Drug 2','FontSize',16); % d2
grid off;
title('Synergy of Potency: CI(Relative LSD)','FontSize',16)
fname_fig = 'syn_potency_CI_LSD_relative';
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);

%% Plot CI(Loewe's): synergy of potency score
figure;
imagesc(d1,d2,CI_Loewes'); 
colorbar
xlabel('Killing rate of Drug 1','FontSize',16); % d1
ylabel('Killing rate of Drug 2','FontSize',16); % d2
grid off;
title('Synergy of Potency: CI(Loewes)','FontSize',16)
fname_fig = 'syn_potency_CI_Loewes';
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
xlabel('Dose (killing rate) of Drug 1','FontSize',16); % d1
ylabel('Dose (killing rate) of Drug 2','FontSize',16); % d2
title('Isobologram','FontSize',16); 
fname_fig = 'isobologram';
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);

%% Pareto front in LSD-HSA multi-synergy criterion space
%CI_LSD_vec = reshape(CI_LSD',1,[]); % Convert matrix to vector data
CI_LSD_vec = reshape(CI_LSD_relEC50',1,[]); % Convert matrix to vector data
CI_HSA_vec = reshape(CI_HSA',1,[]);
fprintf('Pareto Front in LSD-HSA Multi-Synergy Criterion Space:\n'); 
[pareto_HSALSD_LSD_plot,pareto_HSALSD_HSA_plot,onPareto_LSD_HSA] = ...
    find_Pareto(CI_LSD_vec,CI_HSA_vec,CI_LSD_relEC50,CI_HSA,drug1,drug2);
    %find_Pareto(CI_LSD_vec,CI_HSA_vec,CI_LSD,CI_HSA,drug1,drug2);
figure; % plot objective function value at each drug dose tested
plot(CI_LSD_vec,CI_HSA_vec,'ob','MarkerFaceColor','b');  
hold on; 
plot(pareto_HSALSD_LSD_plot,pareto_HSALSD_HSA_plot,'-or','MarkerFaceColor','r'); hold off;
xlabel('Synergy of potency: CI(LSD)','FontSize',16); % objective function 2
ylabel('Synergy of efficacy: CI(HSA)','FontSize',16); % objective function 1
title('Pareto Front (Red)','FontSize',16);
fname_fig = 'pareto_LSD_HSA';
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);

%% Doses on Pareto front in LSD-HSA multi-synergy criterion space
figure; hold on;
imagesc(d1,d2,TGI_combo'); 
colorbar
caxis([0,1]);
xlabel('Killing rate of Drug 1','FontSize',16); % d1
ylabel('Killing rate of Drug 2','FontSize',16); % d2
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
%plot(d1_experiment,d2_experiment,'*r','linewidth',0.5);
hold off; 
xlim([0.9*min(d1),max(d1)+0.1*min(d1)]); 
ylim([0.9*min(d2),max(d2)+0.1*min(d2)]); 
fname_fig = 'doses_onPareto_LSD_HSA';
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);


%% Pareto front in LSD-Bliss multi-synergy criterion space
CI_Bliss_vec = reshape(CI_Bliss',1,[]); % Convert matrix to vector data
fprintf('Pareto Front in LSD-Bliss Multi-Synergy Criterion Space:\n'); 
[pareto_BlissLSD_LSD_plot,pareto_BlissLSD_Bliss_plot,onPareto_LSD_Bliss] = ...
    find_Pareto(CI_LSD_vec,CI_Bliss_vec,CI_LSD_relEC50,CI_Bliss,drug1,drug2);
    %find_Pareto(CI_LSD_vec,CI_Bliss_vec,CI_LSD,CI_Bliss,drug1,drug2);
figure; % plot objective function value at each drug dose tested
plot(CI_LSD_vec,CI_Bliss_vec,'ob','MarkerFaceColor','b');  
hold on; 
plot(pareto_BlissLSD_LSD_plot,pareto_BlissLSD_Bliss_plot,'-or','MarkerFaceColor','r'); hold off;
xlabel('Synergy of potency: CI(LSD)','FontSize',16); % objective function 2
ylabel('Synergy of efficacy: CI(Bliss)','FontSize',16); % objective function 1
title('Pareto Front (Red)','FontSize',16);
fname_fig = 'pareto_LSD_Bliss';
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);

%% Doses on Pareto front in LSD-Bliss multi-synergy criterion space
figure; hold on;
imagesc(d1,d2,TGI_combo'); 
colorbar
caxis([0,1]);
xlabel('Killing rate of Drug 1','FontSize',16); % d1
ylabel('Killing rate of Drug 2','FontSize',16); % d2
grid off;
title('x = on LSD-Bliss Pareto, Color=TGI','FontSize',16);
for i = 1:length(d1)
    for j = 1:length(d2)
        if(onPareto_LSD_Bliss(i,j)>0)
            plot(drug1(i,j),drug2(i,j),'xk','linewidth',2);
        end
    end
end
%plot(d1_experiment,d2_experiment,'*r','linewidth',0.5);
hold off; 
xlim([0.9*min(d1),max(d1)+0.1*min(d1)]); 
ylim([0.9*min(d2),max(d2)+0.1*min(d2)]); 
fname_fig = 'doses_onPareto_LSD_Bliss';
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);

%% Pareto front in Loewes-HSA multi-synergy criterion space
CI_Loewes_vec = reshape(CI_Loewes',1,[]); % Convert matrix to vector data
fprintf('Pareto Front in Loewes-HSA Multi-Synergy Criterion Space:\n'); 
[pareto_LoewesHSA_Loewes_plot,pareto_LoewesHSA_HSA_plot,onPareto_Loewes_HSA] = ...
    find_Pareto(CI_Loewes_vec,CI_HSA_vec,CI_Loewes,CI_HSA,drug1,drug2);
figure; % plot objective function value at each drug dose tested
plot(CI_Loewes_vec,CI_HSA_vec,'ob','MarkerFaceColor','b');  
hold on; 
plot(pareto_LoewesHSA_Loewes_plot,pareto_LoewesHSA_HSA_plot,'-or','MarkerFaceColor','r'); hold off;
xlabel('Synergy of potency: CI(Loewes)','FontSize',16); % objective function 2
ylabel('Synergy of efficacy: CI(HSA)','FontSize',16); % objective function 1
title('Pareto Front (Red)','FontSize',16);
fname_fig = 'pareto_Loewes_HSA';
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);

%% Doses on Pareto front in Loewes-HSA multi-synergy criterion space
figure; hold on;
imagesc(d1,d2,TGI_combo'); 
colorbar
caxis([0,1]);
xlabel('Killing rate of Drug 1','FontSize',16); % d1
ylabel('Killing rate of Drug 2','FontSize',16); % d2
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
%plot(d1_experiment,d2_experiment,'*r','linewidth',0.5);
hold off; 
xlim([0.9*min(d1),max(d1)+0.1*min(d1)]); 
ylim([0.9*min(d2),max(d2)+0.1*min(d2)]); 
fname_fig = 'doses_onPareto_Loewes_HSA';
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);


%% Pareto front in Loewes-Bliss multi-synergy criterion space
fprintf('Pareto Front in Loewes-Bliss Multi-Synergy Criterion Space:\n'); 
[pareto_LoewesBliss_Loewes_plot,pareto_LoewesBliss_Bliss_plot,onPareto_Loewes_Bliss] = ...
    find_Pareto(CI_Loewes_vec,CI_Bliss_vec,CI_Loewes,CI_Bliss,drug1,drug2);
figure; % plot objective function value at each drug dose tested
plot(CI_Loewes_vec,CI_Bliss_vec,'ob','MarkerFaceColor','b');  
hold on; 
plot(pareto_LoewesBliss_Loewes_plot,pareto_LoewesBliss_Bliss_plot,'-or','MarkerFaceColor','r'); hold off;
xlabel('Synergy of potency: CI(Loewes)','FontSize',16); % objective function 2
ylabel('Synergy of efficacy: CI(Bliss)','FontSize',16); % objective function 1
title('Pareto Front (Red)','FontSize',16);
fname_fig = 'pareto_Loewes_Bliss';
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);

%% Doses on Pareto front in Loewes-Bliss multi-synergy criterion space
figure; hold on;
imagesc(d1,d2,TGI_combo'); 
colorbar
caxis([0,1]);
xlabel('Killing rate of Drug 1','FontSize',16); % d1
ylabel('Killing rate of Drug 2','FontSize',16); % d2
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
%plot(d1_experiment,d2_experiment,'*r','linewidth',0.5);
hold off; 
xlim([0.9*min(d1),max(d1)+0.1*min(d1)]); 
ylim([0.9*min(d2),max(d2)+0.1*min(d2)]); 
fname_fig = 'doses_onPareto_Loewes_Bliss';
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
xlabel('Killing rate of Drug 1','FontSize',16); % d1
ylabel('Killing rate of Drug 2','FontSize',16); % d2
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
%plot(d1_experiment,d2_experiment,'*r','linewidth',0.5);
hold off; 
fname_fig = 'doses_onPareto_any';
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

toc
diary off
save output.mat TGI_drug1 TGI_drug2 TGI_combo CI_HSA CI_Bliss closest_D1 ...
    closest_D2 CI_LSD CI_Loewes drug1 drug2 cum_dose min_D pareto_HSALSD_LSD_plot ...
    pareto_HSALSD_HSA_plot pareto_BlissLSD_LSD_plot pareto_BlissLSD_Bliss_plot ...
    pareto_LoewesHSA_Loewes_plot pareto_LoewesHSA_HSA_plot ...
    pareto_LoewesBliss_Loewes_plot pareto_LoewesBliss_Bliss_plot onPareto ...
    onPareto_LSD_HSA onPareto_LSD_Bliss onPareto_Loewes_HSA onPareto_Loewes_Bliss ...
    CI_LSD_relEC50 min_D_relEC50



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions needed to run script                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Find dose of specified monotherapy closest to target TGI
function [closest_D, closest_D_TGI] = closest_drug(d,d_num,TGI_target)
    global r M y0 tf sol_control
    dist_closest_D = 10^10; closest_D = -1; closest_D_TGI = -1; 
    for k = 1:length(d)
        if d_num==1
            sol_d = ode23s(@(t,x) logistic_drug(t,x,r,M,d(k),0),[0,tf],y0);
        elseif d_num==2
            sol_d = ode23s(@(t,x) logistic_drug(t,x,r,M,0,d(k)),[0,tf],y0);
        else
            fprintf('Can only use d_num = 1 or 2, but set d_num = %d\n',d_num);
            STOP
        end
        TGI_d = (sol_control.y(end)-sol_d.y(end))/sol_control.y(end);
        dist_d = abs(TGI_d-TGI_target); 
        if(dist_d<dist_closest_D)
            closest_D = d(k); 
            closest_D_TGI = TGI_d; 
            dist_closest_D = dist_d; 
        end
    end
end

%% Find points on Pareto front
function [par_potency_plot,par_efficacy_plot,pareto] = find_Pareto(CI_potency_vec,...
    CI_efficacy_vec,CI_pot,CI_eff,drug1,drug2)
    global onPareto

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
        fprintf('\tPareto point #%d with CI_pot = %f and CI_eff = %f\n',k,pareto_pot(k),pareto_eff(k)); 
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


%% Toy system of DEs
function xp = logistic_drug(t,x,r,K,d1,d2)
 xp = r*x*(1-x/K)-d1*x-d2*x^2;
end

