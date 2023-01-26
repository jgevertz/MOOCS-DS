clearvars; close all; clc;
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
    path = 'Pembro_Bevacizumab/ForCluster_Schedules/CellLine1_H1299/'; 
    %path = 'ForCluster/Test/'; 
elseif(cell_line == 2)
    path = 'Pembro_Bevacizumab/ForCluster_Schedules/CellLine2_A549/'; 
else
    fprintf('Can only enter 1 or 2 for cell line - exiting\n'); 
    stop 
end
s = [path 'Output_1_1/output1_1.mat']; 
load(s); 

cmap = copper(length(spacing)*length(spacing));
count = 1; 
figure; hold on;
for i = 1:length(spacing)
    for j = 1:length(spacing)
        if ((i==1) && (j==1))
            % Do nothing
        else
            clear pareto_HSALSD_LSD_plot pareto_HSALSD_HSA_plot; 
            s = [path 'Output_' num2str(spacing(i)) '_' num2str(spacing(j))...
                '/output' num2str(spacing(i)) '_' num2str(spacing(j)) '.mat'];
            load(s);
        end
        if j == 1
            h1 = plot(pareto_HSALSD_LSD_plot,pareto_HSALSD_HSA_plot,...
                'LineWidth',2,'Color', cmap(count,:),'DisplayName',...
                sprintf('%d/%d',spacing(i),spacing(j)));
        elseif j == 2
            h1 = plot(pareto_HSALSD_LSD_plot,pareto_HSALSD_HSA_plot,...
                '--','LineWidth',2,'Color', cmap(count,:),'DisplayName',...
                sprintf('%d/%d',spacing(i),spacing(j)));            
        elseif j == 3
            h1 = plot(pareto_HSALSD_LSD_plot,pareto_HSALSD_HSA_plot,...
                '-.','LineWidth',2,'Color', cmap(count,:),'DisplayName',...
                sprintf('%d/%d',spacing(i),spacing(j)));
        elseif j == 4
            h1 = plot(pareto_HSALSD_LSD_plot,pareto_HSALSD_HSA_plot,...
                ':','LineWidth',2,'Color', cmap(count,:),'DisplayName',...
                sprintf('%d/%d',spacing(i),spacing(j)));
        else 
            h1 = plot(pareto_HSALSD_LSD_plot,pareto_HSALSD_HSA_plot,...
                '+-','LineWidth',2,'Color', cmap(count,:),'DisplayName',...
                sprintf('%d/%d',spacing(i),spacing(j)));
        end
        count = count + 1;
    end
end
hold off;
legend('Location','NorthEast','FontSize',14,'NumColumns',3)
xlabel('Synergy of potency: CI(LSD)','FontSize',16); % objective function 2
ylabel('Synergy of efficacy: CI(HSA)','FontSize',16); % objective function 1

s = [path 'Output_1_2/output1_2.mat']; 
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
%plot(d1_experiment,d2_experiment,'*r','linewidth',0.5);
hold off; 
xlim([0.9*min(d1),max(d1)+0.1*min(d1)]); 
ylim([0.9*min(d2),max(d2)+0.1*min(d2)]); 