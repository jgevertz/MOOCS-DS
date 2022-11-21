%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% (10/5/22: JG) Pareto front across protocols (days between doses)      %
% - Modification of the main code which gives drugs Q3D so the analysis %
%   can be done for any protocol (that is, any separation of drug1 and  %
%   drug2).                                                             %
% - In the ForCluster folder, created an m-file for each protocol       %
%   to analyze. Here, that allows for dose of each drug to be separated %
%   by 1 to 5 days. So there are 25 different protocols to consider.    %
% - The entire multi-objective synergy is performed by running each     %
%   code. It's easiest to run each protocol case on the cluster. The    %
%   output is saved to a dynamically named folder. If m is the (fixed)  %
%   spacing between doses of drug1 and n is the (fixed) spacing between %
%   doses of drug2, all output gets saved in the folder Output_m_n.     %
% - Must make sure spacing vector represents the actual protocols ran.  %
%   Other than this, everything is automatic                            %
% - In dosing space, the code displays the TGI of the experimental      %
%   protocol, and labels all points that appear on the Pareto front of  %
%   multiple protocols. 0 values aren't indicated. The max number of    %
%   protocol Pareto fronts a dose can fall on is:                       %
%     length(spacing)*length(spacing)                                   %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set/read in data
clear all; close all; clc;
prompt = "What cell line do you want to analyze? Enter 1 for H1299, 2 for A549: ";
cell_line = input(prompt);
spacing = 1:1:5; % Test case: spacing = 1:2:3; 
onPareto_Loewes_Bliss_all = cell(length(spacing)*length(spacing)); 
onPareto_Loewes_HSA_all = cell(length(spacing)*length(spacing));  
onPareto_LSD_Bliss_all = cell(length(spacing)*length(spacing)); 
onPareto_LSD_HSA_all = cell(length(spacing)*length(spacing)); 
d1_experiment = 10; % actual dose of drug 1
d2_experiment = 1; % actual dose of drug 2
font = 'Arial';
if(cell_line == 1)
    path = 'ForCluster/CellLine1_H1299/'; 
    %path = 'ForCluster/Test/'; 
elseif(cell_line == 2)
    path = 'ForCluster/CellLine2_A549/'; 
else
    fprintf('Can only enter 1 or 2 for cell line - exiting\n'); 
    stop 
end
s = [path 'Output_1_1/output1_1.mat']; 
load(s); 
num_pts = length(drug1);
min_d1 = min(drug1,[],'all');
max_d1 = max(drug1,[],'all');
d1 = linspace(min_d1,max_d1,num_pts);
min_d2 = min(drug2,[],'all');
max_d2 = max(drug2,[],'all');
d2 = linspace(min_d2,max_d2,num_pts);
step_d1 = 0.5*(max(d1)-min(d1))/num_pts; 
step_d2 = 0.5*(max(d2)-min(d2))/num_pts; 

count = 1; 
for i = 1:length(spacing)
    for j = 1:length(spacing)
        if ((i==1) && (j==1))
            % Do nothing
        else
            clear onPareto_Loewes_Bliss onPareto_Loewes_HSA ...
                  onPareto_LSD_Bliss onPareto_LSD_HSA; 
            s = [path 'Output_' num2str(spacing(i)) '_' num2str(spacing(j))...
                '/output' num2str(spacing(i)) '_' num2str(spacing(j)) '.mat'];
            load(s);
        end
        onPareto_Loewes_Bliss_all{count} = onPareto_Loewes_Bliss; 
        onPareto_Loewes_HSA_all{count} = onPareto_Loewes_HSA;  
        onPareto_LSD_Bliss_all{count} = onPareto_LSD_Bliss; 
        onPareto_LSD_HSA_all{count} = onPareto_LSD_HSA; 
        count = count+1; 
    end
end

clear TGI_combo
s = [path 'Output_3_3/output3_3.mat']; 
load(s); % This way plots show TGI_combo from experimental protocol
onPareto_Loewes_Bliss_acrossProtocols = zeros(length(d1),length(d2));
onPareto_Loewes_HSA_acrossProtocols = zeros(length(d1),length(d2));
onPareto_LSD_Bliss_acrossProtocols = zeros(length(d1),length(d2));
onPareto_LSD_HSA_acrossProtocols = zeros(length(d1),length(d2));
for k = 1:length(onPareto_Loewes_Bliss_all)
    for i = 1:length(d1)
        for j = 1:length(d2)
            if(onPareto_Loewes_Bliss_all{k}(i,j)>0)
                onPareto_Loewes_Bliss_acrossProtocols(i,j) = ...
                    onPareto_Loewes_Bliss_acrossProtocols(i,j) + 1; 
                % fprintf('i = %d, j = %d on Pareto front %d\n',i,j,k); 
            end
            if(onPareto_Loewes_HSA_all{k}(i,j)>0)
                onPareto_Loewes_HSA_acrossProtocols(i,j) = ...
                    onPareto_Loewes_HSA_acrossProtocols(i,j) + 1; 
                % fprintf('i = %d, j = %d on Pareto front %d\n',i,j,k); 
            end
            if(onPareto_LSD_Bliss_all{k}(i,j)>0)
                onPareto_LSD_Bliss_acrossProtocols(i,j) = ...
                    onPareto_LSD_Bliss_acrossProtocols(i,j) + 1; 
                % fprintf('i = %d, j = %d on Pareto front %d\n',i,j,k); 
            end
            if(onPareto_LSD_HSA_all{k}(i,j)>0)
                onPareto_LSD_HSA_acrossProtocols(i,j) = ...
                    onPareto_LSD_HSA_acrossProtocols(i,j) + 1; 
                % fprintf('i = %d, j = %d on Pareto front %d\n',i,j,k); 
            end
        end
    end
end
%% Doses on Pareto front in Loewes-Bliss multi-synergy criterion space
fprintf('Finding Pareto Front in Loewes-Bliss Multi-Synergy Criterion Space\n'); 
figure; hold on;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.3, 0.65]); % [left bottom width height]
imagesc(d1,d2,TGI_combo'); 
colorbar
caxis([0,1]);
grid off;
xlabel('Dose pembrolizumab (mg/kg)','FontSize',16); % d1
ylabel('Dose bevacizumab (mg/kg)','FontSize',16); % d2
grid off;
for i = 1:length(d1)
    for j = 1:length(d2)
        if(onPareto_Loewes_Bliss_acrossProtocols(i,j)>0)
            m = num2str(onPareto_Loewes_Bliss_acrossProtocols(i,j)); 
            text(drug1(i,j),drug2(i,j),m,'fontname',font,'color','black','fontweight','bold');
        end
    end
end
plot(d1_experiment,d2_experiment,'*r','linewidth',0.5);
hold off; 
xlim([min(d1)-step_d1,max(d1)+step_d1]); 
ylim([min(d2)-step_d2,max(d2)+step_d2]); 
title('Number of Protocols for which Dose is on Loewe-Bliss Pareto Front','FontSize',14);
subtitle('Color Gives TGI for Experimental Protocol (3 day intervals)','FontSize',14);
fname_fig = [path 'doses_onPareto_Loewes_Bliss'];
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);

%% Doses on Pareto front in Loewes-HSA multi-synergy criterion space
fprintf('Finding Pareto Front in Loewes-HSA Multi-Synergy Criterion Space\n'); 
figure; hold on;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.3, 0.65]); % [left bottom width height]
imagesc(d1,d2,TGI_combo'); 
colorbar
caxis([0,1]);
grid off;
xlabel('Dose pembrolizumab (mg/kg)','FontSize',16); % d1
ylabel('Dose bevacizumab (mg/kg)','FontSize',16); % d2
grid off;
for i = 1:length(d1)
    for j = 1:length(d2)
        if(onPareto_Loewes_HSA_acrossProtocols(i,j)>0)
            m = num2str(onPareto_Loewes_HSA_acrossProtocols(i,j)); 
            text(drug1(i,j),drug2(i,j),m,'fontname',font,'color','black','fontweight','bold');
        end
    end
end
plot(d1_experiment,d2_experiment,'*r','linewidth',0.5);
hold off; 
xlim([min(d1)-step_d1,max(d1)+step_d1]); 
ylim([min(d2)-step_d2,max(d2)+step_d2]); 
title('Number of Protocols for which Dose is on Loewe-HSA Pareto Front','FontSize',14);
subtitle('Color Gives TGI for Experimental Protocol (3 day intervals)','FontSize',14);
fname_fig = [path 'doses_onPareto_Loewes_HSA'];
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);

%% Doses on Pareto front in LSD-Bliss multi-synergy criterion space
fprintf('Finding Pareto Front in LSD-Bliss Multi-Synergy Criterion Space\n'); 
figure; hold on;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.3, 0.65]); % [left bottom width height]
imagesc(d1,d2,TGI_combo'); 
colorbar
caxis([0,1]);
grid off;
xlabel('Dose pembrolizumab (mg/kg)','FontSize',16); % d1
ylabel('Dose bevacizumab (mg/kg)','FontSize',16); % d2
grid off;
for i = 1:length(d1)
    for j = 1:length(d2)
        if(onPareto_LSD_Bliss_acrossProtocols(i,j)>0)
            m = num2str(onPareto_LSD_Bliss_acrossProtocols(i,j)); 
            text(drug1(i,j),drug2(i,j),m,'fontname',font,'color','black','fontweight','bold');
        end
    end
end
plot(d1_experiment,d2_experiment,'*r','linewidth',0.5);
hold off; 
xlim([min(d1)-step_d1,max(d1)+step_d1]); 
ylim([min(d2)-step_d2,max(d2)+step_d2]); 
title('Number of Protocols for which Dose is on LSD-Bliss Pareto Front','FontSize',14);
subtitle('Color Gives TGI for Experimental Protocol (3 day intervals)','FontSize',14);
fname_fig = [path 'doses_onPareto_LSD_Bliss'];
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);

%% Doses on Pareto front in Loewes-HSA multi-synergy criterion space
fprintf('Finding Pareto Front in LSD-HSA Multi-Synergy Criterion Space\n'); 
figure; hold on;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.3, 0.65]); % [left bottom width height]
imagesc(d1,d2,TGI_combo'); 
colorbar
caxis([0,1]);
grid off;
xlabel('Dose pembrolizumab (mg/kg)','FontSize',16); % d1
ylabel('Dose bevacizumab (mg/kg)','FontSize',16); % d2
grid off;
for i = 1:length(d1)
    for j = 1:length(d2)
        if(onPareto_LSD_HSA_acrossProtocols(i,j)>0)
            m = num2str(onPareto_LSD_HSA_acrossProtocols(i,j)); 
            text(drug1(i,j),drug2(i,j),m,'fontname',font,'color','black','fontweight','bold');
        end
    end
end
plot(d1_experiment,d2_experiment,'*r','linewidth',0.5);
hold off; 
xlim([min(d1)-step_d1,max(d1)+step_d1]); 
ylim([min(d2)-step_d2,max(d2)+step_d2]); 
title('Number of Protocols for which Dose is on LSD-HSA Pareto Front','FontSize',14);
subtitle('Color Gives TGI for Experimental Protocol (3 day intervals)','FontSize',14);
fname_fig = [path 'doses_onPareto_LSD_HSA'];
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);


