%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% Reads in output of parameter_uncertainty_3D MOOCS_DS.m (saved in      %
% output.mat) and processes the data. This generates a large number of  %
% figures, as follows:                                                  %
% - near_pareto_per_param_region_individual: Creates one of these for   %
%   each of the Npixels subregions of parameter space. For              %
%   parameterizations in that voxel, displays the probability of a      %
%   combination dose being (near) Pareto optimal in each multi-         %
%   objective synergy space. Doses that do not meet the TGI threshold   %
%   are blacked out.                                                    %
% - near_pareto_per_param_region_all: Creates one of these for each of  %
%   Npixels subregions of parameter space. For parameterizations in     %
%   that voxel, displays the probability of a combination dose being    %
%   Pareto optimal across all multi-objective synergy spaces (the       %
%   Pareto optimality score). Doses that do not meet the TGI threshold  %
%   are blacked out.                                                    %
% - near_pareto_per_param_region_all_TGI: Shows the same probability as %
%   near_pareto_per_param_region_all figures, but shows the probability %
%   on top of a heatmap that gives the TGI associated with each         %
%   combination dose.                                                   %
% - near_pareto_per_combo_dose_all_TGIthresh: Creates one of these for  %
%   each set of combination doses considered. Within each figure,       %
%   displays the Pareto optimality score for each parameterization      %
%   that meets the TGI threshold.                                       %
%                                                                       %
% Data is saved to output_processed_data.mat                            %
% Updated: 8/30/2025                                                    %
% Elapsed time is 464.435645 seconds.                                   %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clear all; tic;
path = 'Output_H1299_3D';
mat_file = [path '/output.mat'];
load(mat_file); 

%% Create a matrix for parameterization k that tells us if a dose is near 
%% Pareto optimal in each of the four criterion spaces
near_dose_matrices_Loewe_Bliss = cell(1, NVPs); % for each VP
near_dose_matrices_Loewe_HSA = cell(1, NVPs);
near_dose_matrices_LSD_Bliss = cell(1, NVPs);
near_dose_matrices_LSD_HSA = cell(1, NVPs); 

for i = 1:NVPs
    near_dose_matrices_Loewe_Bliss{i} = zeros(length(d1),length(d2)); 
    %% Collect near Pareto statistics in Loewe-Bliss: Pareto optimal
    fprintf('Pareto optimal in Loewe-Bliss for VP #%d:\n',i);
    for j = 1: size(pareto_dose_Loewe_Bliss{i},1)
        dose1 = pareto_dose_Loewe_Bliss{i}(j,1);
        dose2 = pareto_dose_Loewe_Bliss{i}(j,2);
        idx_dose1 = find(dose1 == d1); 
        idx_dose2 = find(dose2 == d2); 
        fprintf('\tDose1 = %.2f at index = %d and Dose2 = %f at index = %d\n',...
            dose1,idx_dose1,dose2,idx_dose2);
        near_dose_matrices_Loewe_Bliss{i}(idx_dose1,idx_dose2) = 1; 
    end

    %% Collect near Pareto statistics in Loewe-Bliss: near Pareto optimal
    fprintf('NEAR Pareto optimal in Loewe-Bliss for VP #%d:\n',i);
    for j = 1: size(near_pareto_dose_Loewe_Bliss{i},1)
        dose1 = near_pareto_dose_Loewe_Bliss{i}(j,1);
        dose2 = near_pareto_dose_Loewe_Bliss{i}(j,2);
        idx_dose1 = find(dose1 == d1); 
        idx_dose2 = find(dose2 == d2); 
        fprintf('\tDose1 = %.2f at index = %d and Dose2 = %f at index = %d\n',...
            dose1,idx_dose1,dose2,idx_dose2);
        near_dose_matrices_Loewe_Bliss{i}(idx_dose1,idx_dose2) = 1; 
    end

    %% Collect near Pareto statistics in Loewe-HSA: Pareto optimal
    near_dose_matrices_Loewe_HSA{i} = zeros(length(d1),length(d2)); 
    fprintf('Pareto optimal in Loewe-HSA for VP #%d:\n',i);
    for j = 1: size(pareto_dose_Loewe_HSA{i},1)
        dose1 = pareto_dose_Loewe_HSA{i}(j,1);
        dose2 = pareto_dose_Loewe_HSA{i}(j,2);
        idx_dose1 = find(dose1 == d1); 
        idx_dose2 = find(dose2 == d2); 
        fprintf('\tDose1 = %.2f at index = %d and Dose2 = %f at index = %d\n',...
            dose1,idx_dose1,dose2,idx_dose2);
        near_dose_matrices_Loewe_HSA{i}(idx_dose1,idx_dose2) = 1; 
    end

    %% Collect near Pareto statistics in Loewe-HSA: near Pareto optimal
    fprintf('NEAR Pareto optimal in Loewe-HSA for VP #%d:\n',i);
    for j = 1: size(near_pareto_dose_Loewe_HSA{i},1)
        dose1 = near_pareto_dose_Loewe_HSA{i}(j,1);
        dose2 = near_pareto_dose_Loewe_HSA{i}(j,2);
        idx_dose1 = find(dose1 == d1); 
        idx_dose2 = find(dose2 == d2); 
        fprintf('\tDose1 = %.2f at index = %d and Dose2 = %f at index = %d\n',...
            dose1,idx_dose1,dose2,idx_dose2);
        near_dose_matrices_Loewe_HSA{i}(idx_dose1,idx_dose2) = 1; 
    end

    %% Collect near Pareto statistics in LSD-Bliss: Pareto optimal
    near_dose_matrices_LSD_Bliss{i} = zeros(length(d1),length(d2)); 
    fprintf('Pareto optimal in LSD-Bliss for VP #%d:\n',i);
    for j = 1: size(pareto_dose_LSD_Bliss{i},1)
        dose1 = pareto_dose_LSD_Bliss{i}(j,1);
        dose2 = pareto_dose_LSD_Bliss{i}(j,2);
        idx_dose1 = find(dose1 == d1); 
        idx_dose2 = find(dose2 == d2); 
        fprintf('\tDose1 = %.2f at index = %d and Dose2 = %f at index = %d\n',...
            dose1,idx_dose1,dose2,idx_dose2);
        near_dose_matrices_LSD_Bliss{i}(idx_dose1,idx_dose2) = 1; 
    end

    %% Collect near Pareto statistics in LSD-Bliss: near Pareto optimal
    fprintf('NEAR Pareto optimal in LSD-Bliss for VP #%d:\n',i);
    for j = 1: size(near_pareto_dose_LSD_Bliss{i},1)
        dose1 = near_pareto_dose_LSD_Bliss{i}(j,1);
        dose2 = near_pareto_dose_LSD_Bliss{i}(j,2);
        idx_dose1 = find(dose1 == d1); 
        idx_dose2 = find(dose2 == d2); 
        fprintf('\tDose1 = %.2f at index = %d and Dose2 = %f at index = %d\n',...
            dose1,idx_dose1,dose2,idx_dose2);
        near_dose_matrices_LSD_Bliss{i}(idx_dose1,idx_dose2) = 1; 
    end

    %% Collect near Pareto statistics in Loewe-HSA: Pareto optimal
    near_dose_matrices_LSD_HSA{i} = zeros(length(d1),length(d2)); 
    fprintf('Pareto optimal in LSD-HSA for VP #%d:\n',i);
    for j = 1: size(pareto_dose_LSD_HSA{i},1)
        dose1 = pareto_dose_LSD_HSA{i}(j,1);
        dose2 = pareto_dose_LSD_HSA{i}(j,2);
        idx_dose1 = find(dose1 == d1); 
        idx_dose2 = find(dose2 == d2); 
        fprintf('\tDose1 = %.2f at index = %d and Dose2 = %f at index = %d\n',...
            dose1,idx_dose1,dose2,idx_dose2);
        near_dose_matrices_LSD_HSA{i}(idx_dose1,idx_dose2) = 1; 
    end

    %% Collect near Pareto statistics in Loewe-HSA: near Pareto optimal
    fprintf('NEAR Pareto optimal in LSD-HSA for VP #%d:\n',i);
    for j = 1: size(near_pareto_dose_LSD_HSA{i},1)
        dose1 = near_pareto_dose_LSD_HSA{i}(j,1);
        dose2 = near_pareto_dose_LSD_HSA{i}(j,2);
        idx_dose1 = find(dose1 == d1); 
        idx_dose2 = find(dose2 == d2); 
        fprintf('\tDose1 = %.2f at index = %d and Dose2 = %f at index = %d\n',...
            dose1,idx_dose1,dose2,idx_dose2);
        near_dose_matrices_LSD_HSA{i}(idx_dose1,idx_dose2) = 1; 
    end
end

%% Compute number of times combination dose is Pareto optimal across  
%% all parameterizations within a pixel
near_dose_matrices_Loewe_Bliss_pixel = cell(1, Npixels);%{};
near_dose_matrices_Loewe_HSA_pixel = cell(1, Npixels);
near_dose_matrices_LSD_Bliss_pixel = cell(1, Npixels);
near_dose_matrices_LSD_HSA_pixel = cell(1, Npixels);
near_dose_matrices_all_pixel = cell(1, Npixels);
for i = 1:Npixels
    near_dose_matrices_Loewe_Bliss_pixel{i} = zeros(length(d1),length(d2)); 
    near_dose_matrices_Loewe_HSA_pixel{i} = zeros(length(d1),length(d2)); 
    near_dose_matrices_LSD_Bliss_pixel{i} = zeros(length(d1),length(d2)); 
    near_dose_matrices_LSD_HSA_pixel{i} = zeros(length(d1),length(d2)); 
    near_dose_matrices_all_pixel{i} = zeros(length(d1),length(d2)); 
end

% Convert each pixel element to a number using (Npixels = 5)
map_back = zeros(Npixels,Nparams);
for vp = 1:NVPs
    pixel = params_pixel(vp,:); % Location of VP
    vec_pos = (Nregions-1)^2*(pixel(1)-1) + (Nregions-1)*(pixel(2)-1) + ...
        (pixel(3)-1) + 1;
    map_back(vec_pos,:) = [pixel(1) pixel(2) pixel(3)];
    fprintf('VP %d of %d in pixel (%d, %d, %d) corresponding to vector position %d\n',...
        vp,NVPs,pixel(1),pixel(2),pixel(3),vec_pos);
    
    % % Now, create a matrix that sums all the matrices in that vector position
    near_dose_matrices_Loewe_Bliss_pixel{vec_pos} = ...
        near_dose_matrices_Loewe_Bliss_pixel{vec_pos} + near_dose_matrices_Loewe_Bliss{vp};
    near_dose_matrices_Loewe_HSA_pixel{vec_pos} = ...
        near_dose_matrices_Loewe_HSA_pixel{vec_pos} + near_dose_matrices_Loewe_HSA{vp};
    near_dose_matrices_LSD_Bliss_pixel{vec_pos} = ...
        near_dose_matrices_LSD_Bliss_pixel{vec_pos} + near_dose_matrices_LSD_Bliss{vp};
    near_dose_matrices_LSD_HSA_pixel{vec_pos} = ...
        near_dose_matrices_LSD_HSA_pixel{vec_pos} + near_dose_matrices_LSD_HSA{vp};   
end

for i = 1:Npixels
    near_dose_matrices_all_pixel{i} = near_dose_matrices_Loewe_Bliss_pixel{i} + ...
        near_dose_matrices_Loewe_HSA_pixel{i} + ...
        near_dose_matrices_LSD_Bliss_pixel{i} + ...
        near_dose_matrices_LSD_HSA_pixel{i};
end

font = 'Arial';
%% Display results as a function of the parameter range
TGI_pixel_avg = cell(1, Npixels); % for each parameter "voxel"   
near_dose_matrices_all_pixel_TGI = near_dose_matrices_all_pixel;
for k = 1:Npixels
    TGI_pixel_avg{k} = zeros(length(d1), length(d2));   
end
step_d1 = 0.5*(max(d1)-min(d1))/num_pts; 
step_d2 = 0.5*(max(d2)-min(d2))/num_pts; 
for i = 1:Npixels
    p1_start = param_vec(1,map_back(i,1));
    p1_end = param_vec(1,map_back(i,1)+1);
    p2_start = param_vec(2,map_back(i,2));
    p2_end = param_vec(2,map_back(i,2)+1);
    p3_start = param_vec(3,map_back(i,3)); 
    p3_end = param_vec(3,map_back(i,3)+1);

    p1_avg = (p1_start+p1_end)/2;
    p2_avg = (p2_start+p2_end)/2;
    p3_avg = (p3_start+p3_end)/2;
    TGI_pixel_avg{i} = find_TGI_matrix_fixed_params(p1_avg,p2_avg,p3_avg,...
        d1,d2,tmax,ICs); 
    
    for m = 1:length(d1)
        for n = 1:length(d2)
            if TGI_pixel_avg{i}(m,n) < TGI_thresh % put a black x
                near_dose_matrices_Loewe_Bliss_pixel{i}(m,n) = nan;
                near_dose_matrices_Loewe_HSA_pixel{i}(m,n) = nan;
                near_dose_matrices_LSD_Bliss_pixel{i}(m,n) = nan;
                near_dose_matrices_LSD_HSA_pixel{i}(m,n) = nan;
            end
        end
    end

    figure; 
    set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.6, 0.85]);
    for j = 1:4
        subplot(2,2,j)
        if j==1
            hHM=heatmap(d1,d2,near_dose_matrices_Loewe_Bliss_pixel{i}'/Nsamples); 
            hHM.Title ='Loewe-Bliss';
        elseif j==2
            hHM=heatmap(d1,d2,near_dose_matrices_Loewe_HSA_pixel{i}'/Nsamples); 
            hHM.Title ='Loewe-HSA';
        elseif j==3
            hHM=heatmap(d1,d2,near_dose_matrices_LSD_Bliss_pixel{i}'/Nsamples); 
            hHM.Title ='LSD-Bliss';
        else
            hHM=heatmap(d1,d2,near_dose_matrices_LSD_HSA_pixel{i}'/Nsamples); 
            hHM.Title ='LSD-HSA';
        end
        hHM.YDisplayData=flip(hHM.YDisplayData);
        hHM.XLabel = 'Dose pembrolizumab (mg/kg)';
        hHM.YLabel = 'Dose bevacizumab (mg/kg)';
        hHM.CellLabelFormat = '%.2f';
        hHM.XDisplayLabels = compose('%.2f', str2double(hHM.XDisplayLabels));
        hHM.YDisplayLabels = compose('%.2f', str2double(hHM.YDisplayLabels));
        hHM.ColorLimits = [0 1];

        sgtitle(['Near Pareto Optimality Score: k2_d_1 \in [' ...
            num2str(p1_start) ',' num2str(p1_end) '], k2_d_2 \in [' ...
            num2str(p2_start) ',' num2str(p2_end) , '], k_3 \in [' ...
            num2str(p3_start) ',' num2str(p3_end)  ,']']);
    end
    fname_fig = [path '/near_pareto_per_param_region_individual' num2str(i)];
    saveas(gcf,[fname_fig,'.fig'])
    saveas(gcf,[fname_fig,'.png'])

    figure; hold on;
    set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.35, 0.65]);
    imagesc(d1,d2,TGI_pixel_avg{i}'); 
    xlim([min(d1)-step_d1,max(d1)+step_d1]); 
    ylim([min(d2)-step_d2,max(d2)+step_d2]); 
    xticks(d1)
    yticks(d2)
    xtickformat('%.2f')
    ytickformat('%.2f')
    ax = gca;
    ax.FontSize = 14;
    colorbar
    clim([0,1]);
    xlabel('Dose pembrolizumab (mg/kg)','FontSize',16); % d1
    ylabel('Dose bevacizumab (mg/kg)','FontSize',16); % d2
    grid off;

    title(['Near Pareto Optimal Score Across Synergy Spaces:' ...
        newline 'k2_d_1 \in [' num2str(p1_start) ',' num2str(p1_end) ...
        '], k2_d_2 \in [' num2str(p2_start) ',' num2str(p2_end) , ...
        '], k_3 \in [' num2str(p3_start) ',' num2str(p3_end) ,']'],...
        'FontSize',16);
    subtitle('Color Gives TGI of Combination: Red is 90% TGI Contour','FontSize',14);
    [M,c] = contour(d1,d2,TGI_pixel_avg{i}',[0.9 0.9],'LineColor','r'); %,'ShowText','on');
    c.LineWidth = 3;
    formatSpec = '%3.2f';
    
    for m = 1:length(d1)
        for n = 1:length(d2)
           r = num2str(near_dose_matrices_all_pixel{i}(m,n)/(Ncriterion*Nsamples),formatSpec); 
           %if(dose_matrices_all_pixel{i}(m,n) > 0)
               text(d1(m)-0.25*step_d1,d2(n)-0.1*step_d2,r,'fontname',font,...
                   'color','black','fontweight','bold','FontSize',11);
           %end
            if TGI_pixel_avg{i}(m,n) < TGI_thresh % put a black x
                near_dose_matrices_all_pixel_TGI{i}(m,n) = nan;
            end
        end
    end
    hold off; 
    fname_fig = [path '/near_pareto_per_param_region_all_TGI' num2str(i)];
    saveas(gcf,[fname_fig,'.fig'])
    saveas(gcf,[fname_fig,'.png']) 

    figure;
    set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.35, 0.65]);
    hHM=heatmap(d1,d2,near_dose_matrices_all_pixel_TGI{i}'/(Ncriterion*Nsamples)); 
    hHM.YDisplayData=flip(hHM.YDisplayData);
    hHM.XLabel = 'Dose pembrolizumab (mg/kg)';
    hHM.YLabel = 'Dose bevacizumab (mg/kg)';
    hHM.CellLabelFormat = '%.2f';
    hHM.XDisplayLabels = compose('%.2f', str2double(hHM.XDisplayLabels));
    hHM.YDisplayLabels = compose('%.2f', str2double(hHM.YDisplayLabels));
    hHM.ColorLimits = [0 1];
    hHM.Title = ['Near Pareto Optimality Score Across Synergy Spaces:' ...
        newline 'k2_d_1 \in [' num2str(p1_start) ',' num2str(p1_end) ...
        '], k2_d_2 \in [' num2str(p2_start) ',' num2str(p2_end) ...
        '], k_3 \in [' num2str(p3_start) ',' num2str(p3_end) ']'];
    fname_fig = [path '/near_pareto_per_param_region_all' num2str(i)];
    saveas(gcf,[fname_fig,'.fig'])
    saveas(gcf,[fname_fig,'.png'])
end    

%% Instead, display as as a function of dose
numCombos = num_pts*num_pts; 
p1 = param_vec(1,:);
p2 = param_vec(2,:);
p3 = param_vec(3,:);
% Set grid line color and style
gridColor = [0.6, 0.6, 0.6]; % light gray
gridStyle = '--';

near_param_matrices_all = cell(1, numCombos); % for each combination (d1, d2)      
TGI_pixel_avg_dose = cell(1, numCombos); % for each parameter "pixel"   
for k = 1:numCombos
    near_param_matrices_all{k} = zeros(length(p1)-1, length(p2)-1, length(p3)-1);   
    TGI_pixel_avg_dose{k} = zeros(length(p1)-1, length(p2)-1, length(p3)-1); 
end
p1_avg_vec = zeros(1,length(p1)-1);
p2_avg_vec = zeros(1,length(p2)-1);
p3_avg_vec = zeros(1,length(p3)-1);
for i = 1:length(p1)-1
    p1_avg_vec(i) = (p1(i)+p1(i+1))/2;
end
for j = 1:length(p2)-1
    p2_avg_vec(j) = (p2(j)+p2(j+1))/2;
end
for k = 1:length(p3)-1
    p3_avg_vec(k) = (p3(k)+p3(k+1))/2;
end
[P1_avg, P2_avg, P3_avg] = meshgrid(p1_avg_vec, p2_avg_vec, p3_avg_vec); 
step_p1 = 0.5*(param_vec(1,end)-param_vec(1,1))/(Nregions-1);
step_p2 = 0.5*(param_vec(2,end)-param_vec(2,1))/(Nregions-1);
step_p3 = 0.5*(param_vec(3,end)-param_vec(3,1))/(Nregions-1);
Ccount = 0; 
for i = 1:length(d1)
    for j = 1:length(d2)
        Pcount = 0;
        Ccount = Ccount + 1;
        fprintf('Up to combination %d with d1 = %f, d2 = %f\n',Ccount,...
            d1(i),d2(j));
        for m = 1:length(p1)-1
            for n = 1:length(p2)-1
                for r = 1:length(p3)-1
                    Pcount = Pcount + 1;
                    fprintf('\tUp to parameterization %d with min(p1) = %f', ...
                        Pcount,p1(m));
                    fprintf(' min(p2) = %f, min(p3) = %f\n',p2(n),p3(r));
                    A = near_dose_matrices_all_pixel{Pcount}(i,j);
                    fprintf('\t\tCurrently param_matrices_all{%d}(%d,%d) = %d for param #%d',...
                        Ccount,m,n,near_param_matrices_all{Ccount}(m,n),Pcount);
                    fprintf(' and just found A = %d\n',A);
                    near_param_matrices_all{Ccount}(m,n,r) = ...
                         near_param_matrices_all{Ccount}(m,n,r) + A;
                    fprintf('\t\tparam_matrices_all{%d}(%d,%d,%d) = %d\n',...
                        Ccount,m,n,r,near_param_matrices_all{Ccount}(m,n,r));
                end
            end
        end

        TGI_pixel_avg_dose{Ccount} = find_TGI_matrix_fixed_dose(p1,p2,p3,...
            d1(i),d2(j),tmax,ICs);
        tgi_for_thresholding = permute(TGI_pixel_avg_dose{Ccount}, [2 1 3]);

        % output values for color
        c = near_param_matrices_all{Ccount}/(Ncriterion*Nsamples); 
        % meshgrid reorders axes as (j,i,k) for p1(i), p2(j), p3(k). 
        % So, we must permute A before flattening and using in scatter3
        c_perm = permute(c, [2 1 3]);
        P1_before_thresh = P1_avg(:);
        P2_before_thresh = P2_avg(:);
        P3_before_thresh = P3_avg(:);
        c_before_thresh = c_perm(:);
        tgi_before_thresh = tgi_for_thresholding(:);
        idx_thresh = find(tgi_before_thresh>=TGI_thresh);

        P1_after_thresh = P1_before_thresh(idx_thresh);
        P2_after_thresh = P2_before_thresh(idx_thresh);
        P3_after_thresh = P3_before_thresh(idx_thresh);
        c_after_thresh  = c_before_thresh(idx_thresh);

        figure;
        set(gcf, 'Units', 'Normalized','OuterPosition', ...
            [0.05, 0.05, 0.35, 0.65]);
        % Draw lines in x-direction (varying p1), at each (p2, p3)
        for m = 1:length(p2)
            for n = 1:length(p3)
                line([p1(1), p1(end)], [p2(m), p2(m)], [p3(n), p3(n)], ...
                    'Color', gridColor, 'LineStyle', gridStyle);
            end
        end
        % Draw lines in y-direction (varying p2), at each (p1, p3)
        for m = 1:length(p1)
            for n = 1:length(p3)
                line([p1(m), p1(m)], [p2(1), p2(end)], [p3(n), p3(n)], ...
                    'Color', gridColor, 'LineStyle', gridStyle);
            end
        end
        % Draw lines in z-direction (varying p3), at each (p1, p2)
        for m = 1:length(p1)
            for n = 1:length(p2)
                line([p1(m), p1(m)], [p2(n), p2(n)], [p3(1), p3(end)], ...
                    'Color', gridColor, 'LineStyle', gridStyle);
            end
        end
        scatter3(P1_after_thresh,P2_after_thresh,P3_after_thresh,300,...
            c_after_thresh(:),'filled',"square")
        hold off;
        xlim([param_vec(1,1),param_vec(1,end)]);
        ylim([param_vec(2,1),param_vec(2,end)]);
        zlim([param_vec(3,1),param_vec(3,end)]);
        xlabel('k2_d_1','FontSize',14);
        ylabel('k2_d_2','FontSize',14);
        zlabel('k_3','FontSize',14);
        title(['Near Pareto Optimality Score Across Synergy Spaces:' ...
            newline 'Pembro dose = ' num2str(d1(i)) ' mg/kg; Beva dose = ' ... 
            num2str(d2(j)) ' mg/kg'], 'FontSize',16);
        view(3); 
        colorbar
        clim([0 1])
        fname_fig = [path '/near_pareto_per_combo_dose_all_TGIthresh' ...
            num2str(TGI_thresh) , '_' num2str(Ccount)];
        saveas(gcf,[fname_fig,'.fig'])
        saveas(gcf,[fname_fig,'.png'])
    end
end
toc; 

fsave = [path '/output_processed_data.mat']; 
save(fsave,'TGI_pixel_avg','p1_avg_vec','p2_avg_vec','p3_avg_vec',...
    'TGI_pixel_avg_dose','d1','d2','p1','p2','p3','Ncriterion','Nsamples',...
    'near_dose_matrices_Loewe_Bliss_pixel','near_dose_matrices_Loewe_HSA_pixel',...
    'near_dose_matrices_LSD_Bliss_pixel','near_dose_matrices_LSD_HSA_pixel',...
    'near_dose_matrices_all_pixel','near_param_matrices_all');

function TGI_combo = find_TGI_matrix_fixed_params(p1,p2,p3,d1,d2,tmax,ICs)
    options = odeset('RelTol', 1.0e-10); 
    p = set_parameters();    
    params = [p.lam1 p1 p2 p3]; % use avg of k2_d1, k2_d2, k3; lambda fixed  

    %% Everything is relative to control, so set doses to 0 to find control solution
    %% Since parameters are fixed, this only needs to be determined once
    dose1 = 0; dose2 = 0;
    [tDose1,cDose1,tDose2,cDose2,ti,tUnion,p] = set_protocol(p,tmax,...
        dose1,dose2);
    sol_control = solve_model(params,tDose1,cDose1,tDose2,cDose2,...
        tUnion,ti,p,ICs,options);
    tumor_control = extract_volume(sol_control,...
        0,dose1,dose2); % 0 = no visual; 1 = visualize
    fprintf('lambda = %f, k2d1 %f, k2d2 = %f, k3 = %f\n',params(1),...
        params(2),params(3),params(4));
    fprintf('End control tumor size = %f\n',tumor_control(end))
    %% Now loop over dosing space
    TGI_combo = zeros(length(d1),length(d2));
    for i = 1:length(d1)
        for j = 1:length(d2)
            dose1 = d1(i);
            dose2 = d2(j);
            [tDose1,cDose1,tDose2,cDose2,ti,tUnion,p] = set_protocol(p,...
                tmax,dose1,dose2);
            sol_combo = solve_model(params,tDose1,cDose1,tDose2,cDose2,...
                tUnion,ti,p,ICs,options);
            tumor_combo = extract_volume(sol_combo,0,dose1,dose2); % combination
            TGI_combo(i,j) = (tumor_control(end)-tumor_combo(end))/tumor_control(end);
            fprintf('\tDose1 = %f, Dose2 = %f, TGI = %f\n',dose1,dose2,...
                TGI_combo(i,j));
        end
    end
end

function TGI_combo = find_TGI_matrix_fixed_dose(p1,p2,p3,d1,d2,tmax,ICs)
    options = odeset('RelTol', 1.0e-10); 
    p = set_parameters();    

    %% Now loop over parameter space
    TGI_combo = zeros(length(p1)-1,length(p2)-1,length(p3)-1);
    for i = 1:length(p1)-1
        p1_avg = (p1(i)+p1(i+1))/2;
        for j = 1:length(p2)-1
            p2_avg = (p2(j)+p2(j+1))/2;
                for k = 1:length(p3)-1
                    p3_avg = (p3(k)+p3(k+1))/2;
                    % use avg of k2_d1, k2_d2, k3, lambda fixed  
                    params = [p.lam1 p1_avg p2_avg p3_avg]; 

                    % Since params change, must find control for each parameterization
                    [tDose1,cDose1,tDose2,cDose2,ti,tUnion,p] = set_protocol(p,...
                        tmax,0,0);
                    sol_control = solve_model(params,tDose1,cDose1,tDose2,cDose2,...
                        tUnion,ti,p,ICs,options);
                    tumor_control = extract_volume(sol_control,0,0,0); % 0 = no visual    
        
                    [tDose1,cDose1,tDose2,cDose2,ti,tUnion,p] = set_protocol(p,...
                        tmax,d1,d2);
                    sol_combo = solve_model(params,tDose1,cDose1,tDose2,cDose2,...
                        tUnion,ti,p,ICs,options);
                    tumor_combo = extract_volume(sol_combo,0,d1,d2); % combination
                    TGI_combo(i,j,k) = (tumor_control(end)-tumor_combo(end))/tumor_control(end);
                end
        end
    end
end

function p = set_parameters()
    % Basic options
    p.iv = 0; % SC
  
    %% parameters
    %% PK parameters
    % drug 1 pembro, 2 compartment
    p.V1_1 = 59.6; 
    p.V2_1 = 44.45; % new (mL/kg)
    
    if p.iv == 1
        p.k01_1 = 0;
    elseif p.iv == 0
        p.k01_1 = 0.11; % SC
    end
    
    p.Cl1_1 = 0.61*24*2; 
    p.Cl2_1 = 18.75*24; %new (mL/kg/d)
    
    p.k10_1 = p.Cl1_1/p.V1_1;
    p.k12_1 = p.Cl2_1/p.V1_1;
    p.k21_1 = p.Cl2_1/p.V2_1;
    
    %% %%%% drug 2, bevacizumab, 1 compartment
    
    p.V1_2  = 119; % ml/kg
    p.Cl1_2 = 13.6 ; % ml/kd/day
    p.k10_2 = p.Cl1_2/p.V1_2;
    
    if p.iv == 1
        p.k01_2 = 0;
    elseif p.iv == 0
        p.k01_2 = 1.3; % SC
    end
    
    % tumor parameters
    p.lam1 = 0.24;
    p.K=1300;

    
    % drug efficacy parameters
    p.k2_d1 = 0.15; % effect of pembro 10mpk
    p.k2_d2 = 0.09; % effect of beva 1 mpk
    p.k3 = 0.002; %combination factor
    
    % Not used
    p.alpha = 1; p.beta  = 1;
    p.k1 = .000575; %transitional states
end

function [tDose1,cDose1,tDose2,cDose2,ti,tUnion,p] = set_protocol(p,...
    tmax,dose1,dose2)
    %% Protocol that changes over MOOCS-DS
    p.doseD1 = dose1*10^3/p.V1_1;%/p.MW1;% mpk,
    p.doseD2 = dose2*10^3/p.V1_2;%/p.MW2;% mpk,
    
    ivdoseD1 = p.doseD1;
    ivdoseD2 = p.doseD2;

    %% Protocol that is fixed: User must set this
    d1_spacing = 3; % frequency of drug 1 (experiments: given every 3 days)
    d2_spacing = 3; % frequency of drug 2 (experiments: given every 3 days)
    d1_experiment = 10; % actual dose of drug 1
    d2_experiment = 1; % actual dose of drug 2
    
    % Drug 1 experimental protocol - pembro, 10mpk
    p.dosenumD1 = 5;   % number of doses
    p.freq1=1; % day
    p.intvlD1 = d1_spacing*p.freq1; % dosing interval, days
    p.tStart1 = 17; % drug 1, time of start
    d1_on_off = ones(1,p.dosenumD1); % 1 = all drug on by default

    % Drug 2 experimental protocol - bevacizumab, 1mpk
    p.dosenumD2 = 6;   % number of doses
    p.freq2=1; % day
    p.intvlD2 = d2_spacing*p.freq2; % dosing interval, days
    p.tStart2 = 14; % drug 2, time of start
    d2_on_off = ones(1,p.dosenumD2); % 0 means drug off, 1 means drug on

    tDose1 = zeros(1,p.dosenumD1); cDose1 = zeros(1,p.dosenumD1);
    for i=1:p.dosenumD1
        tDose1(i)=p.intvlD1*(i-1)+ p.tStart1;
        cDose1(i)=ivdoseD1*d1_on_off(i);
    end

    tDose2 = zeros(1,p.dosenumD2); cDose2 = zeros(1,p.dosenumD2);
    for i=1:p.dosenumD2
        tDose2(i)=p.intvlD2*(i-1)+ p.tStart2;
        cDose2(i)=ivdoseD2*d2_on_off(i);
    end

    % Time interval is the UNION of both tDose1 and tDose2
    tUnion = union( tDose1, tDose2, 'sorted');
    if find(tUnion==0) == 1 % do nothing, t = 0 is in there already
    else % Otherwise, include t = 0
        tUnion = [0 tUnion];
    end
    ti = [tUnion tmax];
end

function sol = solve_model(VP_params,tDose1,cDose1,tDose2,cDose2,tUnion,...
    ti,p,ICs,options) 
    p.lam1  = VP_params(1);
    p.k2_d1 = VP_params(2);
    p.k2_d2 = VP_params(3);
    p.k3    = VP_params(4);

    sol = cell(length(tUnion));
    for i = 1:length(tUnion)
        % Does that TIME exist in tDose1 array?
        [lia1,locb1] = ismember( tUnion(i), tDose1 );
        if lia1 == 1
            if p.iv == 1 %IV
                ICs(2) = ICs(2) + cDose1(locb1);  % give drug 1: D1_p0 = y(2)
            elseif p.iv == 0 %SC
                ICs(1) = ICs(1) + cDose1(locb1); % give drug 1: D1_iv0 = y(1)
            end
        end

        % Does that TIME exist in tDose2 array?
        [lia2,locb2] = ismember( tUnion(i), tDose2 );
        if lia2 == 1
            if p.iv == 1 %IV
                ICs(7) = ICs(7) + cDose2(locb2); %give drug 2: D2_p0  = y(7);
            elseif p.iv == 0
                ICs(6) = ICs(6) + cDose2(locb2); %give drug 2: D2_iv0 = y(6);
            end
        end

        % Build the tSpan and run the Diff Eq
        tSpan = [ti(i) ti(i+1)];
        sol{i} = ode45(@(t,y) combination_model(t,y,p),tSpan,ICs,options);

        % Set new ICs and repeat 
        for j = 1:length(ICs)
            ICs(j) = sol{i}.y(j,end); 
        end
    end
end

%% Combination model: pembro+beva
function dydt = combination_model(t,y,p) 
    D1_iv = y(1);
    D1_p  = y(2);
    D1_t  = y(3);
    D1T1  = y(4);
    T1    = y(5);
    
    D2_iv = y(6);
    D2_p  = y(7);
    D2_t  = y(8);
    D2T2  = y(9);
    T2    = y(10);
    
    x1    = y(11);
    x2    = y(12);
    x3    = y(13);
    x4    = y(14);
    wt=x1+x2+x3+x4;

    %% Drug 1: pembro, 2 compartment
    dD1_iv=-p.k01_1*D1_iv;
    dD1_p=p.k01_1*D1_iv-p.k10_1*D1_p...
        -p.k12_1*D1_p+p.k21_1*D1_t*p.V2_1/p.V1_1;
    dD1_t=p.k12_1*D1_p*(p.V1_1/p.V2_1) - p.k21_1*D1_t;
    dT1D1 = 1; %free, to have option to add interaction with target
    dT1 = 1; % %free, to have option to add target turnover

    %% Drug 2: bevacizumab, 1 compartment
    dD2_iv=-p.k01_2*D2_iv; %absorption
    dD2_p=p.k01_2*D2_iv-p.k10_2*D2_p; %distribution and clearance
    dD2_t=1; %free, to have option to do 2 compartment
    dT2D2 = 1; %free, to have option to add target turnover
    dT2 = 1; %

    %%%%%%%%% the important equation %%%%%%%
    kill_term = p.k2_d1*D1_p/(D1_p+25) + p.k2_d2*D2_p/(D2_p+10) + ...
        p.k3*(D1_p)*(D2_p)/((D1_p)*(D2_p) + 15);
    
    %% Tumor
    dx1 = p.lam1*x1*(1-wt/p.K) - x1*kill_term;
    %%%%%%%%% the important equation %%%%%%%

    dx2 = 0*x1*kill_term - p.k1*x2; % remove transition to death compartments
    dx3 = p.k1*(x2 -x3);
    dx4 = p.k1*(x3 -x4);

    dydt=[dD1_iv;dD1_p; dD1_t; dT1D1; dT1;... %drug 1 and target 1
            dD2_iv;dD2_p; dD2_t; dT2D2; dT2;... %drug 2 and targe 2
            dx1 ; dx2 ; dx3; dx4]; %
end

function tumor = extract_volume(sol,to_visualize,dose1,dose2)
        time = []; x1 = []; x2 = []; x3 = []; x4 = []; 
        for j = 1:size(sol,2)
            time = [time sol{j}.x];
            x1 = [x1 sol{j}.y(11,:)]; 
            x2 = [x2 sol{j}.y(12,:)]; 
            x3 = [x3 sol{j}.y(13,:)]; 
            x4 = [x4 sol{j}.y(14,:)]; 
        end
        tumor = x1+x2+x3+x4;
        
        if to_visualize == 1
            figure;
            plot(time,tumor); 
            xlabel('Time','FontSize',14)
            ylabel('Tumor','FontSize',14)
            title(['Growth with d1 = ' num2str(dose1) ', d2 = ' ...
                num2str(dose2)],'FontSize',16);
        end
end