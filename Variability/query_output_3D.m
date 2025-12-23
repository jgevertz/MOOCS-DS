%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% Reads in data processed in analyze_data_3D_TGIthresh.m.               %
% Allows user to query the output to find either:                       %
% 1) Patients (parameter values for k2_d1, k2_d2, k3) that are most     %
%    likely to be Pareto optimal for a specific dose.                   %
% 2) Pareto optimal doses for a specific "patient" (within a given      %
%    range for parameter values k2_d1, k2_d2, k3).                      %
% User indirectly determines how many solutions will be displayed by    %
% setting threshold_for_opt to a desired value. All probabilities       %
% within (threshold_for_opt*100)% of optimally Pareto solution are      %
% provided as output. The full heatmap associated with the specified    %
% dose or patient is also displayed.                                    %
%                                                                       %
% Test cases have been checked to confirm output matches data.          %
%                                                                       %
% Updated: 8/30/2025.                                                   %                                    
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all;
format shortG
threshold_for_opt = 0.1;
path_read = 'Output_H1299_3D';
mat_file = [path_read '/output_processed_data.mat'];
load(mat_file);  
prompt = "Query by dose (enter 1) or patient parameters (enter 2): ";
query_type = input(prompt); % Stores value answered to prompt

% Directory to save figures and output
cl = clock; clN = 0;
for ii = 2:4 % Month/day/hour/minute
    clN = floor(100*clN + cl(ii));
end
if query_type == 1
    path = [path_read '/Doses_' , num2str(clN)];
elseif query_type == 2
    path = [path_read '/Params_' , num2str(clN)];
end
if exist(path, 'dir') ~= 7
    mkdir(path)
end
myfile = "patient_output.txt";
filepath = fullfile(path, myfile);
diary(filepath);

if query_type == 1 % Need to select a dose
    fprintf('Drug 1 (pembro) dose options are:\n')
    disp([1:length(d1); d1]')
    prompt2 = "Enter the integer-value in the first column corresponding to the desired drug 1 dose: ";
    dose1 = input(prompt2); % Stores value answered to prompt

    fprintf('\nDrug 2 (beva) dose options are:\n')
    disp([1:length(d2); d2]')
    prompt3 = "Enter the integer-value in the first column corresponding to the desired drug 2 dose: ";
    dose2 = input(prompt3); % Stores value answered to prompt

    % Convert desired dose into correct location in param_matrices_all
    dose_location = zeros(length(d1),length(d2));
    Ccount = 0;
    for i = 1:length(d1)
        for j = 1:length(d2)
            Ccount = Ccount + 1;
            dose_location(i,j) = Ccount;
        end
    end
    % matrix of interest
    A = near_param_matrices_all{dose_location(dose1,dose2)}; 
    M = max(A,[],'all');

    p1_avg_pos = [];
    p2_avg_pos = [];
    p3_avg_pos = [];
    p1_opt = [];
    p2_opt = [];
    p3_opt = [];
    num_opt = [];
    count = 0; 
    for i = 1:length(p1_avg_vec)
        for j = 1:length(p2_avg_vec)
            for k = 1:length(p3_avg_vec)
                norm_dist_to_M = (M-A(i,j,k))/M;
                % Save all within (100*threshold_for_opt)% of optimal
                if norm_dist_to_M <= threshold_for_opt
                    fprintf('i = %d, j = %d, k = %d, normalized dist = %f\n',...
                        i,j,k,norm_dist_to_M);
                    count = count + 1;
                    p1_avg_pos(count) = i;
                    p2_avg_pos(count) = j;
                    p3_avg_pos(count) = k;
                    p1_opt(count) = p1_avg_vec(i);
                    p2_opt(count) = p2_avg_vec(j);
                    p3_opt(count) = p3_avg_vec(k);
                    num_opt(count) = A(i,j,k);
                end
            end
        end
    end
    
    % Now sort them
    [num_opt_sort,idx] = sort(num_opt,'descend');
    p1_opt_sort = p1_opt(idx);
    p2_opt_sort = p2_opt(idx);
    p3_opt_sort = p3_opt(idx);
    p1_avg_pos_sort = p1_avg_pos(idx);
    p2_avg_pos_sort = p2_avg_pos(idx);
    p3_avg_pos_sort = p3_avg_pos(idx);

    fprintf('Top patients for your combination dose: (pembro, beva) = (%f, %f):\n',...
        d1(dose1),d2(dose2));
    for i = 1:length(num_opt_sort)
        p1_avg = (p1(p1_avg_pos_sort(i)) + p1(p1_avg_pos_sort(i)+1))/2;
        p2_avg = (p2(p2_avg_pos_sort(i)) + p2(p2_avg_pos_sort(i)+1))/2;
        p3_avg = (p3(p3_avg_pos_sort(i)) + p3(p3_avg_pos_sort(i)+1))/2;
        fprintf('\t%d. %f <= k2_d1 <= %f (avg = %f),\n',i,...
            p1(p1_avg_pos_sort(i)), p1(p1_avg_pos_sort(i)+1),p1_avg);
        fprintf('\t%f <= k2_d2 <= %f (avg = %f),\n', ...
            p2(p2_avg_pos_sort(i)), p2(p2_avg_pos_sort(i)+1),p2_avg);
        fprintf('\t%f <= k_3 <= %f (avg = %f) with\n',...
            p3(p3_avg_pos_sort(i)), p3(p3_avg_pos_sort(i)+1),p3_avg);
        fprintf('\t\tProb(Pareto optimal) = %f\n', ...
            num_opt_sort(i)/(Ncriterion*Nsamples))
        %fprintf('\t\tTGI(combo) = \n');
        x = p1_avg_pos_sort(i);
        y = p2_avg_pos_sort(i);
        z = p3_avg_pos_sort(i);
        fprintf('\t\tTGI(combo) = %f\n',...
            TGI_pixel_avg_dose{dose_location(dose1,dose2)}(x,y,z));
    end

    % And visualize all
    tgi_thresh = 0.6;
    Ccount = dose_location(dose1,dose2);
    tgi_for_thresholding = permute(TGI_pixel_avg_dose{Ccount}, [2 1 3]);
    c = near_param_matrices_all{Ccount}/...
        (Ncriterion*Nsamples);  % output values for color

    % So, we must permute A before flattening and using in scatter3
    c_perm = permute(c, [2 1 3]);
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
    P1_before_thresh = P1_avg(:);
    P2_before_thresh = P2_avg(:);
    P3_before_thresh = P3_avg(:);
    c_before_thresh = c_perm(:);
    tgi_before_thresh = tgi_for_thresholding(:);
    idx_thresh = find(tgi_before_thresh>=tgi_thresh);

    P1_after_thresh = P1_before_thresh(idx_thresh);
    P2_after_thresh = P2_before_thresh(idx_thresh);
    P3_after_thresh = P3_before_thresh(idx_thresh);
    c_after_thresh  = c_before_thresh(idx_thresh);
    % Set grid line color and style
    gridColor = [0.6, 0.6, 0.6]; % light gray
    gridStyle = '--';

    figure;
    set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.35, 0.65]);
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
    xlim([p1(1),p1(end)]);
    ylim([p2(1),p2(end)]);
    zlim([p3(1),p3(end)]);
    xlabel('k2_d_1','FontSize',14);
    ylabel('k2_d_2','FontSize',14);
    zlabel('k_3','FontSize',14);
    title(['Near Pareto Optimality Score:' ...
        newline 'Pembro dose = ' num2str(d1(dose1),'%.2f') ...
        ' mg/kg; Beva dose = ' num2str(d2(dose2),'%.2f') ' mg/kg'],...
        'FontSize',16);
    view(3); 
    colorbar
    clim([0 1])
    fname_fig = [path '/pareto_per_combo_dose_TGIthresh' ...
        num2str(tgi_thresh) '_d1_' num2str(dose1) '_d2_' num2str(dose2)];
    saveas(gcf,[fname_fig,'.fig'])
    saveas(gcf,[fname_fig,'.png'])


elseif query_type == 2 % Need to select a parameterization
    fprintf('Parameter 1 (k2_d1) range options are in column 2 and 3:\n')
    disp([1:length(p1_avg_vec); p1(1:end-1); p1(2:end)]')
    prompt2 = "Enter the integer-value in the first column corresponding to the desired range for param 1: ";
    param1 = input(prompt2); % Stores value answered to prompt

    fprintf('Parameter 2 (k2_d2) range options are in column 2 and 3:\n')
    disp([1:length(p2_avg_vec); p2(1:end-1); p2(2:end)]')
    prompt3 = "Enter the integer-value in the first column corresponding to the desired range for param 2: ";
    param2 = input(prompt3); % Stores value answered to prompt

    fprintf('Parameter 3 (k_3) range options are in column 2 and 3:\n')
    disp([1:length(p3_avg_vec); p3(1:end-1); p3(2:end)]')
    prompt4 = "Enter the integer-value in the first column corresponding to the desired range for param 3: ";
    param3 = input(prompt4); % Stores value answered to prompt

    % Convert desired parameterization into correct location in 
    % dose_matrices_all_pixel
    param_location = zeros(length(p1_avg_vec),length(p2_avg_vec),...
        length(p3_avg_vec));
    Pcount = 0;
    for i = 1:length(p1_avg_vec)
        for j = 1:length(p2_avg_vec)
            for k = 1:length(p3_avg_vec)
                Pcount = Pcount + 1;
                param_location(i,j,k) = Pcount;
            end
        end
    end
    param_pos = param_location(param1,param2,param3);
    A = near_dose_matrices_all_pixel{param_pos}; % matrix of interest
    M = max(A,[],'all');

    d1_opt = [];
    d2_opt = [];
    d1_opt_idx = [];
    d2_opt_idx = [];
    num_opt = [];
    count = 0; 
    for i = 1:length(d1)
        for j = 1:length(d2)
            norm_dist_to_M = (M-A(i,j))/M;
            % Save all within (100*threshold_for_opt)% of optimal
            if norm_dist_to_M <= threshold_for_opt
                count = count + 1;
                d1_opt_idx(count) = i;
                d2_opt_idx(count) = j;
                d1_opt(count) = d1(i);
                d2_opt(count) = d2(j);
                num_opt(count) = A(i,j);
            end
        end
    end
    
    % Now sort them
    [num_opt_sort,idx] = sort(num_opt,'descend');
    d1_opt_sort = d1_opt(idx);
    d2_opt_sort = d2_opt(idx);
    d1_opt_idx_sort = d1_opt_idx(idx);
    d2_opt_idx_sort = d2_opt_idx(idx);

    fprintf('Top combination doses for patient with params in range:\n');
    fprintf('%f <= k2_d1 <= %f, %f <= k2_d2 <= %f, %f <= k_3 <= %f\n', ...
        p1(param1), p1(param1+1),p2(param2),p2(param2+1),p3(param3),...
        p3(param3+1))
    for i = 1:length(num_opt_sort)
        fprintf('\t%d. Drug 1 dose = %f mg/kg, drug 2 dose = %f mg/kg:\n',...
            i,d1_opt_sort(i),d2_opt_sort(i));
        fprintf('\t\tProb(Pareto optimal) = %f\n',num_opt_sort(i)/...
            (Ncriterion*Nsamples));
        x = d1_opt_idx_sort(i);
        y = d2_opt_idx_sort(i);
        fprintf('\t\tTGI(combo) = %f\n',...
            TGI_pixel_avg{param_location(param1,param2,param3)}(x,y));
    end

    % And visualize all
    num_pts = length(d1);
    step_d1 = 0.5*(max(d1)-min(d1))/num_pts; 
    step_d2 = 0.5*(max(d2)-min(d2))/num_pts; 
    figure; hold on;
    set(gcf, 'Units', 'Normalized','OuterPosition',...
        [0.05, 0.05, 0.4, 0.75]);
    imagesc(d1,d2,TGI_pixel_avg{param_pos}'); 
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
    
    title(['Near Pareto Optimality Score:' ...
        newline 'k2_d_1 \in [' num2str(p1(param1)) ',' num2str(p1(param1+1)) ...
        '], k2_d_2 \in [' num2str(p2(param2)) ',' num2str(p2(param2+1)) , ...
        '], k_3 \in [' num2str(p3(param3)) ',' num2str(p3(param3+1)) ,']'],...
        'FontSize',16);
    subtitle('Color Gives TGI of Combination: Red is 90% TGI Contour',...
        'FontSize',14);
    [M,c] = contour(d1,d2,TGI_pixel_avg{param_pos}',[0.9 0.9],...
        'LineColor','r'); 
    c.LineWidth = 3;
    formatSpec = '%3.2f';
    font = 'Arial';
    for m = 1:length(d1)
        for n = 1:length(d2)
           r = num2str(near_dose_matrices_all_pixel{param_pos}(m,n)/...
               (Ncriterion*Nsamples),formatSpec); 
           text(d1(m)-0.25*step_d1,d2(n)-0.1*step_d2,r,'fontname',font,...
               'color','black','fontweight','bold','FontSize',11);
        end
    end
    hold off; 
    fname_fig = [path '/pareto_k2d1_' num2str(param1) '_k2d2_' ...
        num2str(param2) '_k3_' num2str(param3)];
    saveas(gcf,[fname_fig,'.fig'])
    saveas(gcf,[fname_fig,'.png'])
end

diary off;
