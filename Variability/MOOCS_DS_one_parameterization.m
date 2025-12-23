%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% MOOCS-DS at a single parameterization. The default is the nominal     %
% parameters, though this can be adjusted as desired. For the other     %
% cases shown in the paper, change:                                     %
%   VPs(vp,2) = 0.22; % gives higher sensitivity to pembro              %
%   VPs(vp,3) = 0.18;  % gives higher sensitivity to beva               %
% Some algorithm hyperparameters:                                       %
% - Set display_graphs to 0 if don't want to see all the graphs         %
% - Set num_pts to determine discretization of dosing space. Paper uses %
%   num_pts = 10 and num_pts = 6 is decent for a test run.              %
% - Removes any combination doses for which the TGI < TGI_thresh.       %
%   Pareto optimization only performed on those doses that meet the TGI %
%   threshold. Paper uses TGI_thresh = 0.6.                             %
% - A point is defined as being "near" Pareto optimal provided it is    %
%   within a 2% relative distance from the Pareto front.                %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc; tic;
Ncriterion = 4; % number of multi-objective optimization spaces
NVPs = 1; % just considering one parameterization
TGI_thresh = 0.6; 
CI_ignore = 10^10;

%% Delcare variables for synergy of output and input scores
TGI_drug1 = cell(NVPs,1); 
TGI_drug2 = cell(NVPs,1); 
TGI_combo = cell(NVPs,1); 
CI_HSA = cell(NVPs,1); 
CI_Bliss = cell(NVPs,1); 
closest_D1 = cell(NVPs,1); 
closest_D2 = cell(NVPs,1); 
CI_LSD_relEC50 = cell(NVPs,1); 
drug1 = cell(NVPs,1); 
drug2 = cell(NVPs,1); 
cum_dose_rel_EC50 = cell(NVPs,1); 
min_D_relEC50 = cell(NVPs,1); 
CI_Loewes = cell(NVPs,1);

%% Create cell arrays for each VP and multi-objective space
pareto_points_LSD_HSA = cell(NVPs,1);
pareto_dose_LSD_HSA = cell(NVPs,1);
pareto_points_LSD_Bliss = cell(NVPs,1);
pareto_dose_LSD_Bliss = cell(NVPs,1);
pareto_points_Loewe_HSA = cell(NVPs,1);
pareto_dose_Loewe_HSA = cell(NVPs,1);
pareto_points_Loewe_Bliss = cell(NVPs,1);
pareto_dose_Loewe_Bliss = cell(NVPs,1);
near_pareto_points_LSD_HSA = cell(NVPs,1);
near_pareto_dose_LSD_HSA = cell(NVPs,1);
near_pareto_points_LSD_Bliss = cell(NVPs,1);
near_pareto_dose_LSD_Bliss = cell(NVPs,1);
near_pareto_points_Loewe_HSA = cell(NVPs,1);
near_pareto_dose_Loewe_HSA = cell(NVPs,1);
near_pareto_points_Loewe_Bliss = cell(NVPs,1);
near_pareto_dose_Loewe_Bliss = cell(NVPs,1);

%% Prepare for MOOCS-DS
options = odeset('RelTol', 1.0e-10); 
display_graphs = 1; % 1 will show all plots 
num_pts = 10; % Want 10 for paper, but 6 is a decent test run
d1_experiment = 10; % actual dose of drug 1
d2_experiment = 1; % actual dose of drug 2
[p,ICs,d1,d2,d1_search,d2_search,min_d1,max_d1,min_d2,max_d2] = ...
    set_parameters_ICs_protocol(num_pts); 

Nparams = 3;
params_all = zeros(NVPs,Nparams);  
%% This runs MOOCS-DS at the nominal parameterization
params_all(1) = p.k2_d1; 
params_all(2) = p.k2_d2;
params_all(3) = p.k3; 

tmax = 35; % Stopping time for simulations
path = 'Output_H1299_3D_Nominal';
if exist(path, 'dir') ~= 7
    mkdir(path)
end
sdiary = [path '/PointsOnPareto_nominal'];
diary(sdiary); 

%% Can change scenario as desired. Scenarios tested in paper:
% If want higher sensitivity to pembro
% VPs(vp,2) = 0.22; % Nominal k2_d1 = 0.15
% If want higher sensitivity to beva:
% VPs(vp,3) = 0.18;  % Nominal k2_d2 = 0.09
% Tumor growth rate, lambda, does not change
vec_front = p.lam1*ones(size(params_all,1),1); 
VPs = [vec_front params_all];

%% Run MOOCS-DS for the 1 VP
vp = 1; 
fprintf('Parameterization: lambda = %f, k2_d_1 = %f, k2_d_2 = %f, k_3 = %f\n',...
    VPs(vp,1),VPs(vp,2),VPs(vp,3),VPs(vp,4));
%% Everything is relative to control, so set doses to 0 to find control solution
%% and use same ICs as combination 
dose1 = 0; dose2 = 0;
[tDose1,cDose1,tDose2,cDose2,ti,tUnion,p] = set_protocol(p,tmax,dose1,dose2);
sol_control = solve_model(VPs(vp,:),tDose1,cDose1,tDose2,cDose2,...
    tUnion,ti,p,ICs,options);
tumor_control = extract_volume(sol_control,1,dose1,dose2); % 0 = no visual; 1 = visualize
fprintf('End control tumor size = %f\n',tumor_control(end))
[tDose1,cDose1,tDose2,cDose2,ti,tUnion,p] = set_protocol(p,tmax,10,1);
sol_combo = solve_model(VPs(vp,:),tDose1,cDose1,tDose2,cDose2,...
            tUnion,ti,p,ICs,options);
tumor_combo = extract_volume(sol_combo,2,10,1); % combination

%% Compute EC50 for each drug
% search over EC5-EC60 range, but a high resolution
d1_EC50 = linspace(min_d1,max_d1,num_pts*num_pts); 
[EC50_d1, EC50_d1_TGI] = closest_drug(d1_EC50,1,0.5,p,tmax,...
    tumor_control,ICs,VPs(vp,:),options);
fprintf('\tEC50 for drug 1 = %f as it has TGI = %f\n',EC50_d1, ...
    EC50_d1_TGI);
d2_EC50 = linspace(min_d2,max_d2,num_pts*num_pts); 
% search over EC5-EC60 range, but a high resolution
[EC50_d2, EC50_d2_TGI] = closest_drug(d2_EC50,2,0.5,p,tmax,...
    tumor_control,ICs,VPs(vp,:),options);
fprintf('\tEC50 for drug 2 = %f as it has TGI = %f\n',EC50_d2,...
    EC50_d2_TGI);

for i = 1:length(d1)
    dose1 = d1(i);  dose2 = 0;
    [tDose1,cDose1,tDose2,cDose2,ti,tUnion,p] = set_protocol(p,...
        tmax,dose1,dose2);
    sol_drug1 = solve_model(VPs(vp,:),tDose1,cDose1,tDose2,cDose2,...
        tUnion,ti,p,ICs,options);
    tumor_drug1 = extract_volume(sol_drug1,0,dose1,dose2); % drug1 mono
    
    for j = 1:length(d2)
        dose1 = 0;  dose2 = d2(j);
        [tDose1,cDose1,tDose2,cDose2,ti,tUnion,p] = set_protocol(p,...
            tmax,dose1,dose2);
        sol_drug2 = solve_model(VPs(vp,:),tDose1,cDose1,tDose2,cDose2,...
            tUnion,ti,p,ICs,options);
        tumor_drug2 = extract_volume(sol_drug2,0,dose1,dose2); % drug2 mono

        dose1 = d1(i);
        [tDose1,cDose1,tDose2,cDose2,ti,tUnion,p] = set_protocol(p,...
            tmax,dose1,dose2);
        sol_combo = solve_model(VPs(vp,:),tDose1,cDose1,tDose2,cDose2,...
            tUnion,ti,p,ICs,options);
        tumor_combo = extract_volume(sol_combo,0,dose1,dose2); % combination
        
        %% Synergy of efficacy (E = TGI = percent cell kill relative to control)
        TGI_drug1{vp}(i,j) = (tumor_control(end)-tumor_drug1(end))/...
            tumor_control(end);
        TGI_drug2{vp}(i,j) = (tumor_control(end)-tumor_drug2(end))/...
            tumor_control(end);
        TGI_combo{vp}(i,j) = (tumor_control(end)-tumor_combo(end))/...
            tumor_control(end);
        fprintf('\tDose(drug1) = %f, Dose(drug2) = %f has TGI(combo) = %f\n',...
            dose1,dose2,TGI_combo{vp}(i,j)); 

        if TGI_combo{vp}(i,j) >= TGI_thresh
            % HSA
            CI_HSA{vp}(i,j) = max(TGI_drug1{vp}(i,j),...
                TGI_drug2{vp}(i,j))/TGI_combo{vp}(i,j);
            % Bliss
            bliss = TGI_drug1{vp}(i,j) + TGI_drug2{vp}(i,j) - ...
                (TGI_drug1{vp}(i,j)*TGI_drug2{vp}(i,j));
            CI_Bliss{vp}(i,j) = bliss/TGI_combo{vp}(i,j);
                   
            %% Synergy of potency
            drug1{vp}(i,j) = d1(i); 
            drug2{vp}(i,j) = d2(j); 
            cum_dose_rel_EC50{vp}(i,j) = d1(i)/EC50_d1+d2(j)/EC50_d2; 
            
            % Find dose of drug 1 that gives TGI_combo
            [closest_D1{vp}(i,j), TGI_closest_D1] = ...
                closest_drug(d1_search,1,TGI_combo{vp}(i,j),...
                p,tmax,tumor_control,ICs,VPs(vp,:),options);  
            fprintf('\t\tTarget TGI = %f: Dose %f of drug 1 has TGI %f\n',...
                TGI_combo{vp}(i,j),closest_D1{vp}(i,j),TGI_closest_D1);
            
            % Find dose of drug 2 that gives TGI_combo
            [closest_D2{vp}(i,j), TGI_closest_D2] = ...
                closest_drug(d2_search,2,TGI_combo{vp}(i,j),...
                p,tmax,tumor_control,ICs,VPs(vp,:),options);
            fprintf('\t\tTarget TGI = %f: Dose %f of drug 2 has TGI %f\n',...
                TGI_combo{vp}(i,j),closest_D2{vp}(i,j),TGI_closest_D2); 
              
            % Loewe
            CI_Loewes{vp}(i,j) = (d1(i)/closest_D1{vp}(i,j)) + ...
                (d2(j)/closest_D2{vp}(i,j));  
            % Relative LSD
            min_D_relEC50{vp}(i,j)= min(closest_D1{vp}(i,j)/EC50_d1, ...
                closest_D2{vp}(i,j)/EC50_d2);
            CI_LSD_relEC50{vp}(i,j) = cum_dose_rel_EC50{vp}(i,j)/...
                min_D_relEC50{vp}(i,j);
        else % ignore, does not meet TGI threshold
            drug1{vp}(i,j) = d1(i); 
            drug2{vp}(i,j) = d2(j); 
            CI_HSA{vp}(i,j) = nan;
            CI_Bliss{vp}(i,j) = nan;
            CI_Loewes{vp}(i,j) = nan;
            CI_LSD_relEC50{vp}(i,j) = nan;
            fprintf('\t\tTGI < %f, ignore dose\n',TGI_thresh);
        end
    end
end

TGI_combo_plot = TGI_combo{vp};
CI_HSA_plot = CI_HSA{vp};
CI_Bliss_plot = CI_Bliss{vp};
cum_dose_rel_EC50_plot = cum_dose_rel_EC50{vp};
CI_LSD_relEC50_plot = CI_LSD_relEC50{vp};
CI_Loewes_plot = CI_Loewes{vp};
if display_graphs == 1
    TGI_drug1_plot = TGI_drug1{vp};
    TGI_drug2_plot = TGI_drug2{vp};
    make_plots_TGI_CI(TGI_drug1_plot,TGI_drug2_plot,TGI_combo_plot,...
        CI_HSA_plot,CI_Bliss_plot,cum_dose_rel_EC50_plot,...
        CI_LSD_relEC50_plot,CI_Loewes_plot,d1,d2,vp);
end

%% Pareto optimization in four multi-objective synergy spaces  
drug1_vec = reshape(drug1{vp}',1,[]);
drug2_vec = reshape(drug2{vp}',1,[]); 
CI_LSD_vec = reshape(CI_LSD_relEC50_plot',1,[]); % Matrix to vector
CI_HSA_vec = reshape(CI_HSA_plot',1,[]);
CI_Bliss_vec = reshape(CI_Bliss_plot',1,[]); 
CI_Loewes_vec = reshape(CI_Loewes_plot',1,[]); 
for obj=1:Ncriterion
    if obj == 1 % Pareto front in LSD-HSA criterion space  
        % Combine vectors into single matrix
        points = [CI_LSD_vec' CI_HSA_vec' drug1_vec' drug2_vec']; 
    elseif obj == 2 % Pareto front in LSD-Bliss criterion space
        points = [CI_LSD_vec' CI_Bliss_vec' drug1_vec' drug2_vec'];
    elseif obj == 3 % Pareto front in Loewe-HSA criterion space
        points = [CI_Loewes_vec' CI_HSA_vec' drug1_vec' drug2_vec'];
    else % Pareto front in Loewe-Bliss criterion space
        points = [CI_Loewes_vec' CI_Bliss_vec' drug1_vec' drug2_vec'];
    end
    sorted_points = sortrows(points, 1); % Sort the points by x-coordinate
    %find Pareto points and corresponding drug doses
    [pareto_points, pareto_dose]=find_pareto(sorted_points); 
    [near_pareto_points, near_pareto_dose] = find_near_pareto(points,...
        pareto_points);
 
    if obj == 1
        pareto_points_LSD_HSA{vp} = pareto_points;
        pareto_dose_LSD_HSA{vp} = pareto_dose;
        near_pareto_points_LSD_HSA{vp} = near_pareto_points;
        near_pareto_dose_LSD_HSA{vp} = near_pareto_dose;
        fprintf('\tPareto Front in LSD-HSA Multi-Synergy Criterion Space:\n');
    elseif obj == 2
        pareto_points_LSD_Bliss{vp} = pareto_points;
        pareto_dose_LSD_Bliss{vp}= pareto_dose;
        near_pareto_points_LSD_Bliss{vp} = near_pareto_points;
        near_pareto_dose_LSD_Bliss{vp} = near_pareto_dose;
        fprintf('Pareto Front in LSD-Bliss Multi-Synergy Criterion Space:\n');
    elseif obj == 3
        pareto_points_Loewe_HSA{vp} = pareto_points;
        pareto_dose_Loewe_HSA{vp} = pareto_dose;
        near_pareto_points_Loewe_HSA{vp} = near_pareto_points;
        near_pareto_dose_Loewe_HSA{vp} = near_pareto_dose;
        fprintf('Pareto Front in Loewe-HSA Multi-Synergy Criterion Space:\n');
    else
        pareto_points_Loewe_Bliss{vp} = pareto_points;
        pareto_dose_Loewe_Bliss{vp} = pareto_dose;
        near_pareto_points_Loewe_Bliss{vp} = near_pareto_points;
        near_pareto_dose_Loewe_Bliss{vp} = near_pareto_dose;
        fprintf('Pareto Front in Loewe-Bliss Multi-Synergy Criterion Space:\n');
    end

    % Displays Pareto Optimal points found in the loop, along with corresponding drug doses
    for k = 1:size(pareto_points, 1)
        fprintf('\t\tPareto point #%d with CI_pot = %.6f and CI_eff = %.6f\n',...
            k, pareto_points(k,1), pareto_points(k,2));
        fprintf('\t\t\tDose: Drug 1 = %.6f, Drug 2 = %.6f\n', ...
            pareto_dose(k,1), pareto_dose(k,2));
    end

    % Display Near Pareto Optimal points
    for k = 1:size(near_pareto_points, 1)
        fprintf('\t\tNear Pareto point #%d with CI_pot = %.6f and CI_eff = %.6f\n', ...
            k, near_pareto_points(k,1), near_pareto_points(k,2));
        fprintf('\t\tDose: Drug 1 = %.6f, Drug 2 = %.6f\n', ...
            near_pareto_dose(k,1), near_pareto_dose(k,2));
    end

    if display_graphs == 1
        make_plots_Pareto(obj,points,pareto_points,pareto_dose,...
            TGI_combo_plot,d1,d2,drug1{vp},drug2{vp},d1_experiment,...
            d2_experiment,vp,near_pareto_points,near_pareto_dose,path); 
    end
end

%% Now count number of times each dose is on or near a Pareto front for this VP
font = 'Arial';
count_onPareto = cell(NVPs,1);
count_nearPareto = cell(NVPs,1);
step_d1 = 0.5*(max(d1)-min(d1))/num_pts; 
step_d2 = 0.5*(max(d2)-min(d2))/num_pts; 
for vp = 1:NVPs
    count_onPareto{vp} = zeros(length(d1),length(d2)); 
    count_nearPareto{vp} = zeros(length(d1),length(d2)); 
    for i = 1:length(d1)
        for j = 1:length(d2)

            %% LSD-HSA
            % These are Pareto optimal (and thus near Pareto as well)
            for k = 1:size(pareto_dose_LSD_HSA{vp}, 1)
                if abs(drug1{vp}(i, j) - pareto_dose_LSD_HSA{vp}(k, 1)) < 1e-6 && ...
                        abs(drug2{vp}(i, j) - pareto_dose_LSD_HSA{vp}(k, 2)) < 1e-6
                    count_onPareto{vp}(i, j) = count_onPareto{vp}(i, j) + 1;
                    count_nearPareto{vp}(i, j) = count_nearPareto{vp}(i, j) + 1;
                end
            end
            % These are only near Pareto optimal
            for k = 1:size(near_pareto_dose_LSD_HSA{vp}, 1)
                if abs(drug1{vp}(i, j) - near_pareto_dose_LSD_HSA{vp}(k, 1)) < 1e-6 && ...
                        abs(drug2{vp}(i, j) - near_pareto_dose_LSD_HSA{vp}(k, 2)) < 1e-6
                    count_nearPareto{vp}(i, j) = count_nearPareto{vp}(i, j) + 1;
                end
            end 

            %% LSD- Bliss
            % These are Pareto optimal (and thus near Pareto as well)
            for k = 1:size(pareto_dose_LSD_Bliss{vp}, 1)
                if abs(drug1{vp}(i, j) - pareto_dose_LSD_Bliss{vp}(k, 1)) < 1e-6 && ...
                        abs(drug2{vp}(i, j) - pareto_dose_LSD_Bliss{vp}(k, 2)) < 1e-6
                    count_onPareto{vp}(i, j) = count_onPareto{vp}(i, j) + 1;
                    count_nearPareto{vp}(i, j) = count_nearPareto{vp}(i, j) + 1;
                end
            end
            % These are only near Pareto optimal
            for k = 1:size(near_pareto_dose_LSD_Bliss{vp}, 1)
                if abs(drug1{vp}(i, j) - near_pareto_dose_LSD_Bliss{vp}(k, 1)) < 1e-6 && ...
                        abs(drug2{vp}(i, j) - near_pareto_dose_LSD_Bliss{vp}(k, 2)) < 1e-6
                    count_nearPareto{vp}(i, j) = count_nearPareto{vp}(i, j) + 1;
                end
            end

            %% Loewe-HSA
            % These are Pareto optimal (and thus near Pareto as well)
            for k = 1:size(pareto_dose_Loewe_HSA{vp}, 1)
                if abs(drug1{vp}(i, j) - pareto_dose_Loewe_HSA{vp}(k, 1)) < 1e-6 && ...
                        abs(drug2{vp}(i, j) - pareto_dose_Loewe_HSA{vp}(k, 2)) < 1e-6
                    count_onPareto{vp}(i, j) = count_onPareto{vp}(i, j) + 1;
                    count_nearPareto{vp}(i, j) = count_nearPareto{vp}(i, j) + 1;
                end
            end
            % These are only near Pareto optimal
            for k = 1:size(near_pareto_dose_Loewe_HSA{vp}, 1)
                if abs(drug1{vp}(i, j) - near_pareto_dose_Loewe_HSA{vp}(k, 1)) < 1e-6 && ...
                        abs(drug2{vp}(i, j) - near_pareto_dose_Loewe_HSA{vp}(k, 2)) < 1e-6
                    count_nearPareto{vp}(i, j) = count_nearPareto{vp}(i, j) + 1;
                end
            end

            %% Loewe-Bliss
            % These are Pareto optimal (and thus near Pareto as well)
            for k = 1:size(pareto_dose_Loewe_Bliss{vp}, 1)
                if abs(drug1{vp}(i, j) - pareto_dose_Loewe_Bliss{vp}(k, 1)) < 1e-6 && ...
                        abs(drug2{vp}(i, j) - pareto_dose_Loewe_Bliss{vp}(k, 2)) < 1e-6
                    count_onPareto{vp}(i, j) = count_onPareto{vp}(i, j) + 1;
                    count_nearPareto{vp}(i, j) = count_nearPareto{vp}(i, j) + 1;
                end
            end
            % These are only near Pareto optimal
            for k = 1:size(near_pareto_dose_Loewe_Bliss{vp}, 1)
                % These are Pareto optimal (and thus near Pareto as well)
                if abs(drug1{vp}(i, j) - near_pareto_dose_Loewe_Bliss{vp}(k, 1)) < 1e-6 && ...
                        abs(drug2{vp}(i, j) - near_pareto_dose_Loewe_Bliss{vp}(k, 2)) < 1e-6
                    count_nearPareto{vp}(i, j) = count_nearPareto{vp}(i, j) + 1;
                end
            end
        end
    end

    if display_graphs == 1
        TGI_combo_plot = TGI_combo{vp};

        %% On Pareto
        figure; hold on;
        imagesc(d1,d2,TGI_combo_plot'); 
        xlim([min(d1)-step_d1,max(d1)+step_d1]); 
        ylim([min(d2)-step_d2,max(d2)+step_d2]); 
        colorbar
        clim([0,1]);
        xlabel('Dose pembrolizumab (mg/kg)','FontSize',16); % d1
        ylabel('Dose bevacizumab (mg/kg)','FontSize',16); % d2
        grid off;
        title(['VP #' num2str(vp) ': Number of Times Dose is Pareto Optimal'],...
            'FontSize',16);
        subtitle('Color Gives TGI of Combination','FontSize',14);
        for i = 1:length(d1)
            for j = 1:length(d2)
               r = num2str(count_onPareto{vp}(i,j));
               if(count_onPareto{vp}(i,j) > 0)
                   text(drug1{vp}(i,j),drug2{vp}(i,j),r,'fontname',font,...
                       'color','black','fontweight','bold');
               end

               if TGI_combo{vp}(i,j) < TGI_thresh % put a black x
                   plot(drug1{vp}(i,j),drug2{vp}(i,j),'xk',...
                       'linewidth',2,'markersize',10);
               end
            end
        end
        plot(d1_experiment,d2_experiment,'*r','linewidth',1,...
            'MarkerSize',10);
        hold off; 

        %% Near Pareto
        figure; hold on;
        imagesc(d1,d2,TGI_combo_plot'); 
        xlim([min(d1)-step_d1,max(d1)+step_d1]); 
        ylim([min(d2)-step_d2,max(d2)+step_d2]); 
        colorbar
        clim([0,1]);
        xlabel('Dose pembrolizumab (mg/kg)','FontSize',16); % d1
        ylabel('Dose bevacizumab (mg/kg)','FontSize',16); % d2
        grid off;
        title(['VP #' num2str(vp) ...
            ': Number of Times Dose is NEAR Pareto Optimal'],...
            'FontSize',16);
        subtitle('Color Gives TGI of Combination','FontSize',14);
        for i = 1:length(d1)
            for j = 1:length(d2)
               r = num2str(count_nearPareto{vp}(i,j));
               if(count_nearPareto{vp}(i,j) > 0)
                   text(drug1{vp}(i,j),drug2{vp}(i,j),r,'fontname',font,...
                       'color','black','fontweight','bold');
               end

               if TGI_combo{vp}(i,j) < TGI_thresh % put a black x
                   fprintf('Should put an x at drug 1 = %f, drug 2 = %f\n',...
                       drug1{vp}(i,j),drug2{vp}(i,j));
                   plot(drug1{vp}(i,j),drug2{vp}(i,j),'xk',...
                       'linewidth',2,'markersize',10);
               end
            end
        end
        plot(d1_experiment,d2_experiment,'*r','linewidth',1,...
            'MarkerSize',10);
        hold off; 
    end
end

%% Average number of times a combination dose is Pareto optimal across VPs
count_onPareto_all = zeros(length(d1),length(d2)); 
count_nearPareto_all = zeros(length(d1),length(d2)); 
for vp = 1:NVPs
    count_onPareto_all = count_onPareto_all + count_onPareto{vp};
    count_nearPareto_all = count_nearPareto_all + count_nearPareto{vp}; 
end
count_onPareto_mean = count_onPareto_all/NVPs;
count_nearPareto_mean = count_nearPareto_all/NVPs;
for i = 1:length(d1)
    for j = 1:length(d2)
        if TGI_combo{vp}(i,j) < TGI_thresh % put a black x
             count_onPareto_mean(i,j) = nan;
             count_nearPareto_mean(i,j) = nan;
        end
    end
end         

figure; 
hHM=heatmap(d1,d2,count_onPareto_mean');
hHM.YDisplayData=flip(hHM.YDisplayData);
hHM.XLabel = 'Dose pembrolizumab (mg/kg)';
hHM.YLabel = 'Dose bevacizumab (mg/kg)';
hHM.Title ='Mean Number of Times Dose is Pareto Optimal Across VPs';
fname_fig = [path '/doses_onPareto'];
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png'])

figure; 
hHM=heatmap(d1,d2,count_nearPareto_mean');
hHM.YDisplayData=flip(hHM.YDisplayData);
hHM.XLabel = 'Dose pembrolizumab (mg/kg)';
hHM.YLabel = 'Dose bevacizumab (mg/kg)';
hHM.Title ='Mean Number of Times Dose is Near Pareto Optimal across VPs';
fname_fig = [path '/doses_nearPareto'];
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png'])

toc
fsave = [path '/output_nominal.mat']; 
save(fsave,'VPs','TGI_drug1','TGI_drug2','TGI_combo','CI_HSA',...
    'CI_Bliss','closest_D1','closest_D2','CI_Loewes','drug1','drug2',...
    'pareto_points_LSD_HSA','pareto_dose_LSD_HSA','pareto_points_LSD_Bliss',...
    'pareto_dose_LSD_Bliss','pareto_points_Loewe_HSA','pareto_dose_Loewe_HSA',...
    'pareto_points_Loewe_Bliss','pareto_dose_Loewe_Bliss','count_onPareto',...
    'count_onPareto_mean','d1','d2','params_all','NVPs',...
    'num_pts','p','tmax','ICs','Ncriterion','near_pareto_points_LSD_HSA',...
    'near_pareto_dose_LSD_HSA','near_pareto_points_LSD_Bliss',...
    'near_pareto_dose_LSD_Bliss','near_pareto_points_Loewe_HSA',...
    'near_pareto_dose_Loewe_HSA','near_pareto_points_Loewe_Bliss',...
    'near_pareto_dose_Loewe_Bliss','count_nearPareto',...
    'count_nearPareto_mean','TGI_thresh');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Functions %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p,ICs,d1,d2,d1_search,d2_search,min_d1,max_d1,min_d2,max_d2] = ...
    set_parameters_ICs_protocol(num_pts)
    % Basic options
    p.iv = 0; % SC

    %% initial conditions
    D1_iv0 = 0; D1_p0  = 0; D1_t0  = 0; % drug 1
    D2_iv0 = 0; D2_p0  = 0;
    x20    = 0; x30    = 0; x40    = 0; %"transitional tumor volume things"
    x10 = 3.5/3; %initial tumor volume

    % these are currently irrelevant, free variables for potential future
    % modifications
    D2_t0  = 0; D2T20  = 0; %drug 2
    D1T10  = 0; %drug 1
    T10    = 0; %concentration of target 1
    T20    = 0; %concentration of target 2
    
    %% parameters
    %% PK parameters
    % drug 1 pembro, 2 compartment
    p.V1_1 = 59.6; 
    p.V2_1 = 44.45; % new (mL/kg)
    p.Cl1_1 = 0.61*24*2; 
    p.Cl2_1 = 18.75*24; %new (mL/kg/d)
    
    if p.iv == 1
        p.k01_1 = 0;
    elseif p.iv == 0
        p.k01_1 = 0.11; % SC
    end
       
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
    p.alpha = 1; p.beta = 1;
    p.k1 = 0.000575;

    ICs = [D1_iv0 D1_p0 D1_t0 D1T10 T10 D2_iv0 D2_p0 D2_t0 D2T20 T20 ...
        x10 x20 x30 x40];

    %% **User must define range of doses to consider. We set the smallest 
    %% dose (min_d1, min_d2) to be the PI5 dose (5% inhibition relative to 
    %% control). We set the largest monotherapy dose (max_d1, max_d2) to be 
    %% the PI70 dose. For this model, the result is that the max efficacy 
    %% of the combination therapy is 95%. 

    % User-defined dosing range to consider for MOOCS-DS: H1299 here
    min_d1 = 0.408; max_d1 = 10.88; % Experiments use 10 mg/kg of pembro (up to EC70)
    d1 = linspace(min_d1,max_d1,num_pts); % EC5 to EC70 for drug 1 monotherapy
    min_d2 = 0.1372; max_d2 = 3.457; % Experiments use 1 mpk bevacizumab (up to EC70) 
    d2 = linspace(min_d2,max_d2,num_pts); % EC5 to EC70 for drug 2 monotherapy
    %% From EC70 range: E(max_d1,max_d2)=0.9495, which is achieved if go 
    %% up to monotherapy D1=120, D2=15
    % Dose range for drug 1 that ensures we can get monotherapy efficacy = E(max_d1,max_d2)
    d1_search = linspace(min(d1),120,30*num_pts); 
    % Dose range for drug 2 that ensures we can get monotherapy efficacy = E(max_d1,max_d2)  
    d2_search = linspace(min(d2),15,30*num_pts);  
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

%% Solve model
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

    dx2 = 0*x1*kill_term - p.k1*x2;
    dx3 = p.k1*(x2 -x3);
    dx4 = p.k1*(x3 -x4);

    dydt=[dD1_iv;dD1_p; dD1_t; dT1D1; dT1;... %drug 1 and target 1
            dD2_iv;dD2_p; dD2_t; dT2D2; dT2;... %drug 2 and targe 2
            dx1 ; dx2 ; dx3; dx4]; %
end

%% Find dose of specified monotherapy closest to target TGI
function [closest_D, closest_D_TGI] = closest_drug(d,d_num,TGI_target,p,...
    tmax,sol_control,ICs,VP_params,options)
    dist_closest_D = 10^10; closest_D = -1; closest_D_TGI = -1; 
    for k = 1:length(d)
        if d_num==1            
            dose1 = d(k);  dose2 = 0;           
        elseif d_num==2
            dose1 = 0; dose2 = d(k);  
        else
            fprintf('Can only use d_num = 1 or 2, but set d_num = %d\n'...
                ,d_num);
            STOP
        end

        [tDose1,cDose1,tDose2,cDose2,ti,tUnion,p] = set_protocol(p,tmax,...
            dose1,dose2);
        sol = solve_model(VP_params,tDose1,cDose1,tDose2,cDose2,...
            tUnion,ti,p,ICs,options);
        
        tumor = extract_volume(sol,0,dose1,dose2); % 1 to visualize
        TGI_d = (sol_control(end)-tumor(end))/sol_control(end);
        dist_d = abs(TGI_d-TGI_target); 
        if(dist_d<dist_closest_D)
            closest_D = d(k); 
            closest_D_TGI = TGI_d; 
            dist_closest_D = dist_d; 
        end
    end
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
            data.H1299.control.days   = [9 12 15 18 21 24 27 30];
            data.H1299.control.volume = [28 83 133 221 400 528 696 1070];
            figure;
            plot(time,tumor); hold on;
            scatter(data.H1299.control.days, data.H1299.control.volume, ...
                50,'o','DisplayName','H1299 Data: Control');  
            hold off;
            xlabel('Time','FontSize',14)
            ylabel('Tumor','FontSize',14)
            title(['Growth with d1 = ' num2str(dose1) ', d2 = ' ...
                num2str(dose2)],'FontSize',16);
            legend('Location', 'northwest', 'FontSize', 12);
        elseif to_visualize == 2
            data.H1299.combo.days     = [9 12 15 18 21 24 27 30];
            data.H1299.combo.volume   = [22 55 71 81 92 108 129 184];
            figure;
            plot(time,tumor); hold on;
            scatter(data.H1299.combo.days,   data.H1299.combo.volume, ...
                50,'o','DisplayName','H1299 Data: Pembro + Beva'); 
            hold off;
            xlabel('Time','FontSize',14)
            ylabel('Tumor','FontSize',14)
            title(['Growth with d1 = ' num2str(dose1) ', d2 = ' ...
                num2str(dose2)],'FontSize',16);
            legend('Location', 'northwest', 'FontSize', 12);
        end
end

%% Find Points on Pareto Front
function [pareto_points, pareto_dose] = find_pareto(sorted_points)
   best_y = Inf;
   % Initialize empty matrices to store Pareto optimal points and drug doses
   pareto_points = []; 
   pareto_dose = [];
   % Iterate over every point in sorted_points matrix
   for i = 1:size(sorted_points, 1) 
        % Check if the current point is Pareto optimal based on y value
       if sorted_points(i, 2) < best_y
           % Save current point to matrices if current point is Pareto optimal
           pareto_points = [pareto_points; sorted_points(i,1:2)]; 
           pareto_dose = [pareto_dose; sorted_points(i,3:4)];
       end
       best_y = min(best_y, sorted_points(i, 2)); % Update best_y
   end
end

function [] = make_plots_TGI_CI(TGI_drug1_plot,TGI_drug2_plot,...
    TGI_combo_plot,CI_HSA_plot,CI_Bliss_plot,cum_dose_rel_EC50_plot,...
    CI_LSD_relEC50_plot,CI_Loewes_plot,d1,d2,vp)

    %% TGI for drug 1
    figure;  hold on;
    imagesc(d1,d2,TGI_drug1_plot'); hold off;
    colorbar
    clim([0,1]);
    xlabel('Dose pembrolizumab (mg/kg)','FontSize',16); % d1
    ylabel('Dose bevacizumab (mg/kg)','FontSize',16); % d2
    grid off;
    title(['VP #' num2str(vp) ': TGI(drug1) Relative to Control'],...
        'FontSize',16)
    
    %% TGI for drug 2
    figure; hold on;
    imagesc(d1,d2,TGI_drug2_plot'); hold off;
    colorbar
    clim([0,1]);
    xlabel('Dose pembrolizumab (mg/kg)','FontSize',16); % d1
    ylabel('Dose bevacizumab (mg/kg)','FontSize',16); % d2
    grid off;
    title(['VP #' num2str(vp) ': TGI(drug2) Relative to Control'],...
        'FontSize',16)
    
    %% Plot TGI for combination therapy
    figure; hold on;
    imagesc(d1,d2,TGI_combo_plot');  hold off
    xlim([0.9*min(d1),max(d1)+0.1*min(d1)]); 
    ylim([0.9*min(d2),max(d2)+0.1*min(d2)]);   
    colorbar
    clim([0,1]);
    xlabel('Dose pembrolizumab (mg/kg)','FontSize',16); % d1
    ylabel('Dose bevacizumab (mg/kg)','FontSize',16); % d2
    grid off;
    title(['VP #' num2str(vp) ': Combination TGI Relative to Control'],...
        'FontSize',16)
    
    %% Plot CI(HSA): synergy of efficacy score
    % Use consistent colorbar across synergy of efficacy heatmaps
    HSA_min = min(CI_HSA_plot,[],'all');
    HSA_max = max(CI_HSA_plot,[],'all');
    Bliss_min = min(CI_Bliss_plot,[],'all'); 
    Bliss_max = max(CI_Bliss_plot,[],'all');
    CI_eff_min = min(HSA_min,Bliss_min);
    CI_eff_max = max(HSA_max,Bliss_max);
    figure; hold on;
    imagesc(d1,d2,CI_HSA_plot'); hold off
    xlim([0.9*min(d1),max(d1)+0.1*min(d1)]); 
    ylim([0.9*min(d2),max(d2)+0.1*min(d2)]);   
    colorbar
    clim([CI_eff_min,CI_eff_max]);
    xlabel('Dose pembrolizumab (mg/kg)','FontSize',16); % d1
    ylabel('Dose bevacizumab (mg/kg)','FontSize',16); % d2
    grid off;
    title(['VP #' num2str(vp) ': Synergy of Efficacy: CI(HSA)'],...
        'FontSize',16)
    
    %% Plot CI(Bliss): synergy of efficacy score
    figure; hold on;
    imagesc(d1,d2,CI_Bliss_plot'); hold off
    xlim([0.9*min(d1),max(d1)+0.1*min(d1)]); 
    ylim([0.9*min(d2),max(d2)+0.1*min(d2)]);  
    colorbar
    clim([CI_eff_min,CI_eff_max]);
    xlabel('Dose pembrolizumab (mg/kg)','FontSize',16); % d1
    ylabel('Dose bevacizumab (mg/kg)','FontSize',16); % d2
    grid off;
    title(['VP #' num2str(vp) ': Synergy of Efficacy: CI(Bliss)'],...
        'FontSize',16)
    
    %% Plot cumulative drug dose relative to EC50
    figure; hold on;
    imagesc(d1,d2,cum_dose_rel_EC50_plot'); hold off
    xlim([0.9*min(d1),max(d1)+0.1*min(d1)]); 
    ylim([0.9*min(d2),max(d2)+0.1*min(d2)]);   
    colorbar
    xlabel('Dose pembrolizumab (mg/kg)','FontSize',16); % d1
    ylabel('Dose bevacizumab (mg/kg)','FontSize',16); % d2
    grid off;
    title(['VP #' num2str(vp) ': Cumulative Drug Dose Relative to EC50'],...
        'FontSize',16)    
    
    %% Plot CI(LSD rel to EC50): synergy of potency score
    % Use consistent colorbar across synergy of potency heatmaps
    LDSrel_min = min(CI_LSD_relEC50_plot,[],'all');
    LDSrel_max = max(CI_LSD_relEC50_plot,[],'all');
    Loewe_min = min(CI_Loewes_plot,[],'all'); 
    Loewe_max = max(CI_Loewes_plot,[],'all');
    CI_pot_min = min(LDSrel_min,Loewe_min);
    CI_pot_max = max(LDSrel_max,Loewe_max);
    figure; hold on;
    imagesc(d1,d2,CI_LSD_relEC50_plot'); hold off
    xlim([0.9*min(d1),max(d1)+0.1*min(d1)]); 
    ylim([0.9*min(d2),max(d2)+0.1*min(d2)]);   
    colorbar
    clim([CI_pot_min,CI_pot_max]);
    xlabel('Dose pembrolizumab (mg/kg)','FontSize',16); % d1
    ylabel('Dose bevacizumab (mg/kg)','FontSize',16); % d2
    grid off;
    title(['VP #' num2str(vp) ': Synergy of Potency: CI(LSD)'],...
        'FontSize',16)
    
    %% Plot CI(Loewe's): synergy of potency score
    figure; hold on;
    imagesc(d1,d2,CI_Loewes_plot'); hold off
    xlim([0.9*min(d1),max(d1)+0.1*min(d1)]); 
    ylim([0.9*min(d2),max(d2)+0.1*min(d2)]);   
    colorbar
    clim([CI_pot_min,CI_pot_max]);
    xlabel('Dose pembrolizumab (mg/kg)','FontSize',16); % d1
    ylabel('Dose bevacizumab (mg/kg)','FontSize',16); % d2
    grid off;
    title(['VP #' num2str(vp) ': Synergy of Potency: CI(Loewe)'],...
        'FontSize',16)
    
    %% Isobolograms
    [X,Y] = meshgrid(d1,d2);
    TGI_target_vec = 0.1:0.1:0.9; 
    figure;
    [M,c]=contourf(X,Y,TGI_combo_plot',TGI_target_vec,'ShowText','on');
    c.LineWidth = 3;
    clabel(M,c,'FontSize',15);
    colorbar
    xlabel('Dose pembrolizumab (mg/kg)','FontSize',16); % d1
    ylabel('Dose bevacizumab (mg/kg)','FontSize',16); % d2
    title(['VP #' num2str(vp) ': Isobologram'],'FontSize',16); 
end

function [] = make_plots_Pareto(obj,points,pareto_points,pareto_dose,...
    TGI_combo_plot,d1,d2,drug1,drug2,d1_experiment,d2_experiment,vp,...
    near_pareto_points,near_pareto_dose,path)

    %% Plot all points in blue and then Pareto optimal points in red
    figure; hold on;
    scatter(points(:,1), points(:,2), 'filled', 'MarkerFaceColor', 'b'); 
    scatter(pareto_points(:,1), pareto_points(:,2), 'filled', ...
        'MarkerFaceColor', 'r'); 
    plot(pareto_points(:,1), pareto_points(:,2), '-r', 'LineWidth', 1);
    if isempty(near_pareto_points) == 0 % not empty, plot
        scatter(near_pareto_points(:,1), near_pareto_points(:,2),60,'s',...
            'filled','MarkerFaceColor', 'r','MarkerEdgeColor','k',...
            'LineWidth',1.5);
    end
    legend ('Non-optimal points', 'Pareto optimal points','Pareto Front',...
        'Near Pareto points','Location', 'Best','FontSize',14);
    title(['VP #' num2str(vp) ': Pareto Front (Red)'],'FontSize',16);
    if obj == 1
        xlabel('Synergy of potency: CI(LSD)','FontSize',16); 
        ylabel('Synergy of efficacy: CI(HSA)','FontSize',16); 
        fname_fig = 'pareto_LSD_HSA';
    elseif obj == 2
        xlabel('Synergy of potency: CI(LSD)','FontSize',16); 
        ylabel('Synergy of efficacy: CI(Bliss)','FontSize',16); 
        fname_fig = 'pareto_LSD_Bliss';
    elseif obj == 3
        xlabel('Synergy of potency: CI(Loewe)','FontSize',16); 
        ylabel('Synergy of efficacy: CI(HSA)','FontSize',16); 
        fname_fig = 'pareto_Loewe_HSA';
    else
        xlabel('Synergy of potency: CI(Loewe)','FontSize',16); 
        ylabel('Synergy of efficacy: CI(Bliss)','FontSize',16); 
        fname_fig = 'pareto_LSD_HSA';
    end
    a = get(gca,'YTickLabel'); % or use XTickLabel, does same thing 
    set(gca,'YTickLabel',a,'FontSize',13)
    fname_fig_pareto = [path '/' fname_fig];
    saveas(gcf,[fname_fig_pareto,'.fig'])
    saveas(gcf,[fname_fig_pareto,'.png'])

    %% Doses on Pareto front 
    figure; hold on;
    imagesc(d1,d2,TGI_combo_plot'); 
    colorbar
    clim([0,1]);
    xlabel('Dose pembrolizumab (mg/kg)','FontSize',16); % d1
    ylabel('Dose bevacizumab (mg/kg)','FontSize',16); % d2
    grid off;
    for k = 1:size(pareto_dose, 1)
        plot(pareto_dose(k,1), pareto_dose(k,2), 'xk', ...
            'linewidth', 2); 
    end
    for k = 1:size(near_pareto_dose, 1)
        plot(near_pareto_dose(k,1), near_pareto_dose(k,2), 'ok', ...
            'linewidth', 2); 
    end

    plot(d1_experiment,d2_experiment,'*r','linewidth',0.5);
    hold off; 
    xlim([0.9*min(d1),max(d1)+0.1*min(d1)]); 
    ylim([0.9*min(d2),max(d2)+0.1*min(d2)]); 
    if obj == 1
        title(['VP #' num2str(vp) ': x = on LSD-HSA Pareto, Color=TGI'],...
            'FontSize',16);
        fname_fig = 'doses_onPareto_LSD_HSA';
    elseif obj == 2
        title(['VP #' num2str(vp) ': x = on LSD-Bliss Pareto, Color=TGI'],...
            'FontSize',16);
        fname_fig = 'doses_onPareto_LSD_Bliss';
    elseif obj == 3
        title(['VP #' num2str(vp) ': x = on Loewe-HSA Pareto, Color=TGI'],...
            'FontSize',16);
        fname_fig = 'doses_onPareto_Loewe_HSA';
    else
        title(['VP #' num2str(vp) ': x = on Loewe-Bliss Pareto, Color=TGI'],...
            'FontSize',16);
        fname_fig = 'doses_onPareto_LSD_HSA';
    end
    fname_fig_doses = [path '/' fname_fig];
    saveas(gcf,[fname_fig_doses,'.fig'])
    saveas(gcf,[fname_fig_doses,'.png'])
end

%% Find Points Near Pareto Front
function [near_pareto_points, near_pareto_dose] = find_near_pareto(points,...
    pareto_points)
    near_pareto_points = [];
    near_pareto_dose = [];
    
    if size(pareto_points,1) > 1 % more than one point on Pareto front
        min_diff = min(diff(pareto_points(:, 1)));
        x_interp = min(pareto_points(:, 1)):min_diff/10: max(pareto_points(:, 1));
        y_interp = interp1(pareto_points(:, 1), pareto_points(:, 2), ...
            x_interp, 'linear');
        
        % Identify Pareto points
        pareto_idx = ismember(points(:, 1:2), pareto_points, 'rows'); 
        
        non_pareto_points = points(~pareto_idx, 1:2); % Non-Pareto points (objective values)
        non_pareto_doses = points(~pareto_idx, 3:4); % Doses corresponding to non-Pareto points
    
        % Determine "size" of criterion space
        pareto_dist_to_origin = zeros(length(pareto_points),1);
        for i = 1:length(pareto_points)
            pareto_dist_to_origin(i) = pdist2(pareto_points(i,:),[0 0]);
        end
        [dist_min, I_min] = min(pareto_dist_to_origin,[],"all"); % closest to origin
        
        non_pareto_dist_to_origin = zeros(length(non_pareto_points),1);
        for i = 1:length(non_pareto_points)
            non_pareto_dist_to_origin(i) = pdist2(non_pareto_points(i,:),[0 0]);
        end
        [dist_max, I_max] = max(non_pareto_dist_to_origin,[],"all"); % further from origin
        
        distance_extreme = pdist2(pareto_points(I_min,:),non_pareto_points(I_max,:));
        dist_threshold = 0.02*distance_extreme;
        fprintf('Extent of region = %f so using threshold = %f\n',...
            distance_extreme,dist_threshold);
    
        distances = pdist2(non_pareto_points, [x_interp', y_interp']); % Calculate distances
        min_distances = min(distances, [], 2); % Minimum distances for non-Pareto points

        % Loop through each non-Pareto point and check near Pareto front
        for i = 1:size(non_pareto_points, 1)
             % Check if the non-Pareto point is within the threshold
            if min_distances(i) < dist_threshold 
                % Save the near-Pareto point
                near_pareto_points = [near_pareto_points; non_pareto_points(i, :)];
                % Save corresponding doses for near-Pareto point
                near_pareto_dose = [near_pareto_dose; non_pareto_doses(i, :)]; 
            end
        end
    end
end

