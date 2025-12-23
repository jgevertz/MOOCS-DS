%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% Finds parameter values for drug kill terms that cover a +-50% range   %
% from nominal parameterization for cell line H1299.                    %
% - k2_d1 range is determined solving model for drug 1 (pembro only)    %
%   and finding the k2_d1 values so that V(tmax) is +-50% compared to   %
%   V(tmax) using the nominal value of k2_d1.                           %
% - k2_d2 range is determined solving model for drug 2 (beva only)      %
%   and finding the k2_d2 values so that V(tmax) is +-50% compared to   %
%   V(tmax) using the nominal value of k2_d2.                           %
% - k3 range is determined solving model for both drugs and finding the %
%   k3 values so that V(tmax) is +-50% compared to V(tmax) using the    %
%   nominal value of k3. Note this looks ONLY at the effect of k3 by    %
%   zeroing out the monotherapy parameters (k2_d1 = k2_d2 = 0).         %
%                                                                       %
% Updated: 12/23/2025                                                   %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

%% Set some variables
options = odeset('RelTol', 1.0e-10);
cell_line_name = 'H1299';
p_default = initialize_parameters();
y0 = zeros(14,1); % initial condition vector. Only need to change
                  % x1(0) = y0(11) to get the correct starting tumor size
y0_guess = 3.5;
tmax = 35;
% resolution for finding best parameter for target tumor size - 1000
num_pts = 1000; 
data = load_experimental_data(); % experimental data
params_to_vary = {'k2_d1', 'k2_d2', 'k3'};
param_ranges = zeros(3,2);
fprintf('\n--- Sensitivity Analysis: %s ---\n', cell_line_name);

%% Effect of drug 1 without drug 2
j = 1;
pname = params_to_vary{j};
% Tumor volume time course at nominal parameterization
p_d1 = initialize_parameters();
p_d1.dose2 = 0; % no drug 2
p_d1.doseD2 = p_d1.dose2 * 1e3 / p_d1.V1_2;
y0_d1 = y0; 
y0_d1(11) = y0_guess/2; % back-calculated tumor size at t = 0
[tDose1,cDose1,tDose2,cDose2,ti,tUnion] = set_protocol(p_d1,tmax);
sol_d1 = solve_model(tDose1,cDose1,tDose2,cDose2,tUnion,ti,p_d1,y0_d1,...
    options);
[time_d1, wt_d1] = extract_volume(sol_d1); 
wt_d1_end = wt_d1(end); % end volume at nomimal parameters

% Now vary k2_d1 until target end volume is met
p0 = p_d1.(pname);     
for_incr = floor(log10(p0));
p_test_d1 = logspace(for_incr-3,for_incr+3,num_pts);    
final_wt_d1 = zeros(size(p_test_d1));
dist_target_up_d1 = zeros(size(p_test_d1));
dist_target_down_d1 = zeros(size(p_test_d1));
target_up = wt_d1_end * 0.5; % smaller tumor with larger parameter
target_down = wt_d1_end * 1.5; % larger tumor with smaller parameter
fprintf('Nominal %s = %g → Tumor = %g\n', pname, p0, wt_d1_end);
fprintf('Searching over the range: [%g, %g]\n',p_test_d1(1),p_test_d1(end));
p_d1_test = p_d1;
for k = 1:length(p_test_d1)
    p_d1_test.(pname) = p_test_d1(k);
    sol = solve_model(tDose1,cDose1,tDose2,cDose2,tUnion,ti,p_d1_test,...
        y0_d1,options);
    [time, wt] = extract_volume(sol);
    final_wt_d1(k) = wt(end); 
    % fprintf('  Param = %g → Tumor = %g (↓ target = %g, ↑ target = %g)\n',...
    %     p_test_d1(k), final_wt_d1(k), target_down, target_up);
    dist_target_up_d1(k) = abs(final_wt_d1(k)-target_up);
    dist_target_down_d1(k) = abs(final_wt_d1(k)-target_down);
end
[min_up_d1, idx_up_d1] = min(dist_target_up_d1);
[min_down_d1, idx_down_d1] = min(dist_target_down_d1);
param_ranges(j,1) = p_test_d1(idx_down_d1);
param_ranges(j,2) = p_test_d1(idx_up_d1);
fprintf('\tClosest to ↓ target = %g is %s = %g with tumor = %g\n',...
    target_down,pname,p_test_d1(idx_down_d1),final_wt_d1(idx_down_d1));
if param_ranges(j,1) == p_test_d1(1)
    fprintf('\t\t*Cannot achieve target: parameter will take on the lower bound of search range*\n');
else
    ratio = p_default.k2_d1/param_ranges(j,1); 
end
fprintf('\tClosest to ↑ target = %g is %s = %g with tumor = %g\n',...
    target_up,pname,p_test_d1(idx_up_d1),final_wt_d1(idx_up_d1));
if param_ranges(j,2) == p_test_d1(end)
    fprintf('\t\t*Cannot achieve target: parameter will take on the upper bound of search range*\n');
end


%% Effect of drug 2 without drug 1
j = j + 1; 
pname = params_to_vary{j};
% Tumor volume time course at nominal parameterization
p_d2 = initialize_parameters();
p_d2.dose1 = 0; % no drug 1
p_d2.doseD1 = p_d2.dose1 * 1e3 / p_d2.V1_1;
y0_d2 = y0; 
y0_d2(11) = y0_guess;
[tDose1,cDose1,tDose2,cDose2,ti,tUnion] = set_protocol(p_d2,tmax);
sol_d2 = solve_model(tDose1,cDose1,tDose2,cDose2,tUnion,ti,p_d2,...
    y0_d2,options);
[time_d2, wt_d2] = extract_volume(sol_d2);
wt_d2_end = wt_d2(end);% end volume at nomimal parameters

% Now vary k2_d2 until target end volume is met
p0 = p_d2.(pname);     
for_incr = floor(log10(p0));
p_test_d2 = logspace(for_incr-3,for_incr+3,num_pts);    
final_wt_d2 = zeros(size(p_test_d2));
dist_target_up_d2 = zeros(size(p_test_d2));
dist_target_down_d2 = zeros(size(p_test_d2));
target_up = wt_d2_end * 0.5; % smaller tumor with larger parameter
target_down = wt_d2_end * 1.5; % larger tumor with smaller parameter
fprintf('Nominal %s = %g → Tumor = %g\n', pname, p0, wt_d2_end);
fprintf('Searching over the range: [%g, %g]\n',p_test_d2(1),p_test_d2(end));
p_d2_test = p_d2;
for k = 1:length(p_test_d2)
    p_d2_test.(pname) = p_test_d2(k);
    sol = solve_model(tDose1,cDose1,tDose2,cDose2,tUnion,ti,p_d2_test,...
        y0_d2,options);
    [time, wt] = extract_volume(sol);
    final_wt_d2(k) = wt(end); 
    % fprintf('  Param = %g → Tumor = %g (↓ target = %g, ↑ target = %g)\n',...
    %     p_test_d2(k), final_wt_d2(k), target_down, target_up);
    dist_target_up_d2(k) = abs(final_wt_d2(k)-target_up);
    dist_target_down_d2(k) = abs(final_wt_d2(k)-target_down);
end
[min_up_d2, idx_up_d2] = min(dist_target_up_d2);
[min_down_d2, idx_down_d2] = min(dist_target_down_d2);
param_ranges(j,1) = p_test_d2(idx_down_d2);
param_ranges(j,2) = p_test_d2(idx_up_d2);
fprintf('\tClosest to ↓ target = %g is %s = %g with tumor = %g\n',...
    target_down,pname,p_test_d2(idx_down_d2),final_wt_d2(idx_down_d2));
if param_ranges(j,1) == p_test_d2(1)
    fprintf('\t\t*Cannot achieve target: parameter will take on the lower bound of search range*\n');
    % manually set to agreed-upon value of 1e-4 for H1299
    % lower bound double distance from nominal than that for k2_d1
    param_ranges(j,1) = p_default.k2_d2/(ratio*2); 
    fprintf('\t\t\tSetting to agreed upon value of %f\n',param_ranges(j,1));
end
fprintf('\tClosest to ↑ target = %g is %s = %g with tumor = %g\n',...
    target_up,pname,p_test_d2(idx_up_d2),final_wt_d2(idx_up_d2));
if param_ranges(j,2) == p_test_d2(end)
    fprintf('\t\t*Cannot achieve target: parameter will take on the upper bound of search range*\n');
end


%% Effect of combination without monotherapy effects
j = j + 1; 
pname = params_to_vary{j};
% Tumor volume time course at nominal parameterization
p_combo = initialize_parameters();
p_combo.k2_d1 = 0;
p_combo.k2_d2 = 0; 
y0_combo = y0; 
y0_combo(11) = y0_guess/3; 
[tDose1,cDose1,tDose2,cDose2,ti,tUnion] = set_protocol(p_combo,tmax);
sol_combo = solve_model(tDose1,cDose1,tDose2,cDose2,tUnion,ti,...
    p_combo,y0_combo,options);
[time_combo, wt_combo] = extract_volume(sol_combo);
wt_combo_end = wt_combo(end); % end volume at nomimal parameters

% Now vary k3 until target end volume is met
p0 = p_combo.(pname);     
for_incr = floor(log10(p0));
p_test_combo = logspace(for_incr-3,for_incr+3,num_pts);    
%p_test_combo = [0 p_test_combo]; %add 0 for extreme lower bound 
final_wt_combo = zeros(size(p_test_combo));
dist_target_up_combo = zeros(size(p_test_combo));
dist_target_down_combo = zeros(size(p_test_combo));
target_up = wt_combo_end * 0.5; % smaller tumor with larger parameter
target_down = wt_combo_end * 1.5; % larger tumor with smaller parameter
fprintf('Nominal %s = %g → Tumor = %g\n', pname, p0, wt_combo_end);
fprintf('Searching over the range: [%g, %g]\n',p_test_combo(1),p_test_combo(end));
p_combo_test = p_combo;
for k = 1:length(p_test_combo)
    p_combo_test.(pname) = p_test_combo(k);
    sol = solve_model(tDose1,cDose1,tDose2,cDose2,tUnion,ti,...
        p_combo_test,y0_combo,options);
    [time, wt] = extract_volume(sol);
    final_wt_combo(k) = wt(end); 
    % fprintf('  Param = %g → Tumor = %g (↓ target = %g, ↑ target = %g)\n',...
    %     p_test_d2(k), final_wt_d2(k), target_down, target_up);
    dist_target_up_combo(k) = abs(final_wt_combo(k)-target_up);
    dist_target_down_combo(k) = abs(final_wt_combo(k)-target_down);
end
[min_up_combo, idx_up_combo] = min(dist_target_up_combo);
[min_down_combo, idx_down_combo] = min(dist_target_down_combo);
param_ranges(j,1) = p_test_combo(idx_down_combo);
param_ranges(j,2) = p_test_combo(idx_up_combo);
fprintf('\tClosest to ↓ target = %g is %s = %g with tumor = %g\n',...
    target_down,pname,p_test_combo(idx_down_combo),final_wt_combo(idx_down_combo));
if param_ranges(j,1) == p_test_combo(1)
    % lower bound double distance from nominal than that for k2_d1
    param_ranges(j,1) = p_default.k3/(ratio*2);
    fprintf('\t\t*Cannot achieve target*\n');
    fprintf('\t\t\tSetting to agreed upon value of %f\n',param_ranges(j,1));
end
fprintf('\tClosest to ↑ target = %g is %s = %g with tumor = %g\n',...
    target_up,pname,p_test_combo(idx_up_combo),final_wt_combo(idx_up_combo));
if param_ranges(j,2) == p_test_combo(end)
    fprintf('\t\t*Cannot achieve target: parameter will take on the upper bound of search range*\n');
end


%% Plot everything
% ---------- Simulate CONTROL (both doses = 0) ----------
p_control = initialize_parameters();
p_control.dose1 = 0;
p_control.dose2 = 0;
p_control.doseD1 = p_control.dose1 * 1e3 / p_control.V1_1;
p_control.doseD2 = p_control.dose2 * 1e3 / p_control.V1_1;
y0_control = y0; 
y0_control(11) = y0_guess;
[tDose1,cDose1,tDose2,cDose2,ti,tUnion] = set_protocol(p_control,tmax);
sol_control = solve_model(tDose1,cDose1,tDose2,cDose2,tUnion,ti,...
    p_control,y0_control,options);
[time_control, wt_control] = extract_volume(sol_control);

% ---------- Simulate COMBO (default doses) ----------
p_mono_combo = initialize_parameters();
y0_mono_combo = y0; 
y0_mono_combo(11) = y0_guess/3; 
[tDose1,cDose1,tDose2,cDose2,ti,tUnion] = set_protocol(p_mono_combo,tmax);
sol_mono_combo = solve_model(tDose1,cDose1,tDose2,cDose2,tUnion,ti,...
    p_mono_combo,y0_mono_combo,options);
[time_mono_combo, wt_mono_combo] = extract_volume(sol_mono_combo);

cmap = lines(5); 
figure; 
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.3, 0.6]);
hold on;
% Plot model-predicted tumor volume
plot(time_control, wt_control,'LineWidth',2,'DisplayName',...
    'Model: Control');  % control model
plot(time_d1, wt_d1, 'LineWidth',2,'DisplayName',...
    'Model: Pembro Only');  % drug 1 monotherapy
plot(time_d2, wt_d2, 'LineWidth',2,'DisplayName',...
    'Model: Beva Only');  % drug 2 monotherapy
plot(time_combo,wt_combo,'LineWidth', 2,'DisplayName',...
    'Model: Combination Effect Only');  % combo effect only
plot(time_mono_combo,wt_mono_combo,'LineWidth', 2,'DisplayName',...
    'Model: Full Pembro + Beva'); % mono + combo (full model)

% Plot model at nominal parameters compared to data
scatter(data.H1299.control.days, data.H1299.control.volume, 50, ...
    cmap(1,:), 'o','DisplayName','Data: Control');  % control model
scatter(data.H1299.pembro.days, data.H1299.pembro.volume, 50, ...
    cmap(2,:), 's','DisplayName','Data: Pembro Only');  % drug 1 mono
scatter(data.H1299.beva.days, data.H1299.beva.volume, 50, ...
    cmap(3,:), 'd','DisplayName','Data: Beva Only');  % drug 2 monotherapy
scatter(data.H1299.combo.days,   data.H1299.combo.volume, 50,  ...
    cmap(5,:), 'x','DisplayName','Data: Full Pembro + Beva'); % mono+combo 
title('H1299: Nominal Parameterization','fontsize',16);
fname1 = 'fits_H1299';
xlabel('Time (days)','fontsize',14);
ylabel('Tumor Volume (mm^3)','fontsize',14);
legend('Location', 'northwest','fontsize',14);
grid on;
ax = gca;
ax.FontSize = 14;
% saveas(gcf,[fname1,'.fig'])
% saveas(gcf,[fname1,'.png'])

%% Now plot nominal compared to parameter boundaries identified
p_d1_up = p_d1;
p_d1_up.(params_to_vary{1}) = p_test_d1(idx_up_d1);
[tDose1,cDose1,tDose2,cDose2,ti,tUnion] = set_protocol(p_d1_up,tmax);
sol_p_d1_up = solve_model(tDose1,cDose1,tDose2,cDose2,tUnion,ti,...
    p_d1_up,y0_d1,options);
[time_d1_up, wt_d1_up] = extract_volume(sol_p_d1_up);

p_d1_down = p_d1;
p_d1_down.(params_to_vary{1}) = p_test_d1(idx_down_d1);
[tDose1,cDose1,tDose2,cDose2,ti,tUnion] = set_protocol(p_d1_down,tmax);
sol_p_d1_down = solve_model(tDose1,cDose1,tDose2,cDose2,tUnion,ti,...
    p_d1_down,y0_d1,options);
[time_d1_down, wt_d1_down] = extract_volume(sol_p_d1_down);

p_d2_up = p_d2;
p_d2_up.(params_to_vary{2}) = p_test_d2(idx_up_d2);
[tDose1,cDose1,tDose2,cDose2,ti,tUnion] = set_protocol(p_d2_up,tmax);
sol_p_d2_up = solve_model(tDose1,cDose1,tDose2,cDose2,tUnion,ti,...
    p_d2_up,y0_d2,options);
[time_d2_up, wt_d2_up] = extract_volume(sol_p_d2_up);

p_d2_down = p_d2;
p_d2_down.(params_to_vary{2}) = p_test_d2(idx_down_d2);
[tDose1,cDose1,tDose2,cDose2,ti,tUnion] = set_protocol(p_d2_down,tmax);
sol_p_d2_down = solve_model(tDose1,cDose1,tDose2,cDose2,tUnion,ti,...
    p_d2_down,y0_d2,options);
[time_d2_down, wt_d2_down] = extract_volume(sol_p_d2_down);

p_combo_up = p_combo;
p_combo_up.(params_to_vary{3}) = p_test_combo(idx_up_combo);
[tDose1,cDose1,tDose2,cDose2,ti,tUnion] = set_protocol(p_combo_up,tmax);
sol_p_combo_up = solve_model(tDose1,cDose1,tDose2,cDose2,tUnion,ti,...
    p_combo_up,y0_combo,options);
[time_combo_up, wt_combo_up] = extract_volume(sol_p_combo_up);

p_combo_down = p_combo;
p_combo_down.(params_to_vary{3}) = p_test_combo(idx_down_combo);
[tDose1,cDose1,tDose2,cDose2,ti,tUnion] = set_protocol(p_combo_down,tmax);
sol_p_combo_down = solve_model(tDose1,cDose1,tDose2,cDose2,tUnion,ti,...
    p_combo_down,y0_combo,options);
[time_combo_down, wt_combo_down] = extract_volume(sol_p_combo_down);

figure; 
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.8, 0.6]);
sgtitle('H1299','fontsize',16,'fontweight','bold');
fname2 = 'param_ranges_H1299';

subplot(1,3,1)  % drug 1 monotherapy
hold on;
plot(time_d1, wt_d1, 'LineWidth',2,'DisplayName',...
    sprintf('k2_d_1 = %0.4g',p_d1.k2_d1)); 
scatter(data.H1299.pembro.days, data.H1299.pembro.volume, 50, ...
    cmap(1,:), 's','DisplayName','Data');  % drug 1 monotherapy
title('Pembro Only','FontSize',14);
plot(time_d1_down, wt_d1_down, 'LineWidth',2,'DisplayName',...
    sprintf('k2_d_1 = %0.4g',p_test_d1(idx_down_d1))); 
plot(time_d1_up, wt_d1_up, 'LineWidth',2,'DisplayName',...
    sprintf('k2_d_1 = %0.4g',p_test_d1(idx_up_d1))); 
hold off;
legend('Location', 'northwest','fontsize',14);
grid on;
ax = gca;
ax.FontSize = 14;
xlim([0,tmax])

subplot(1,3,2) % drug 2 monotherapy
hold on;
plot(time_d2, wt_d2, 'LineWidth',2,'DisplayName',...
    sprintf('k2_d_2 = %0.4g',p_d1.k2_d2)); 
scatter(data.H1299.beva.days, data.H1299.beva.volume, 50, ...
    cmap(1,:), 'd','DisplayName','Data');  % drug 2 monotherapy
title('Beva Only','FontSize',14);
plot(time_d2_down, wt_d2_down, 'LineWidth',2,'DisplayName',...
    sprintf('k2_d_2 = %0.4g',p_test_d2(idx_down_d2))); 
plot(time_d2_up, wt_d2_up, 'LineWidth',2,'DisplayName',...
    sprintf('k2_d_2 = %0.4g',p_test_d2(idx_up_d2))); 
legend('Location', 'northwest','fontsize',14);
grid on;
ax = gca;
ax.FontSize = 14;
xlim([0,tmax])

subplot(1,3,3) % combo effect only (no data for this)
hold on;
plot(time_combo,wt_combo,'LineWidth', 2,'DisplayName',...
    sprintf('k_3 = %0.4g',p_d1.k3)); 
title('Combination Effect Only (k2_d_1 = k2_d_2 = 0)','FontSize',14);
plot(time_combo_down, wt_combo_down, 'LineWidth',2,'DisplayName',...
    sprintf('k_3 = %0.4g',p_test_combo(idx_down_combo))); 
plot(time_combo_up, wt_combo_up, 'LineWidth',2,'DisplayName',...
    sprintf('k_3 = %0.4g',p_test_combo(idx_up_combo))); 
legend('Location', 'northwest','fontsize',14);
grid on;
ax = gca;
ax.FontSize = 14;
xlim([0,tmax])
saveas(gcf,[fname2,'.fig'])
saveas(gcf,[fname2,'.png'])

fsave = [fname2 '.mat']; 
save(fsave,'param_ranges'); % row 1: k2_d1
                            % row 2: k2_d2
                            % row 3: k_3


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p = initialize_parameters()
    % Basic options
    p.iv = 0; % SC

    % Pembro
    p.tStart1 = 17; % drug 1, time of start
    p.dose1 = 10;
    p.dosenumD1 = 5;
    p.freq1 = 1;
    p.intvlD1 = 3 * p.freq1;
    p.V1_1 = 59.6; 
    p.V2_1 = 44.45; % new (mL/kg)
    p.Cl1_1 = 0.61*24*2; 
    p.Cl2_1 = 18.75*24; %new (mL/kg/d)
    p.k01_1 = 0.11;
    p.k10_1 = p.Cl1_1 / p.V1_1;
    p.k12_1 = p.Cl2_1 / p.V1_1;
    p.k21_1 = p.Cl2_1 / p.V2_1;
    p.doseD1 = p.dose1 * 1e3 / p.V1_1;

    % Beva
    p.tStart2 = 14; % drug 2, time of start
    p.dose2 = 1;
    p.dosenumD2 = 6;
    p.freq2 = 1;
    p.intvlD2 = 3 * p.freq2;
    p.V1_2 = 119;
    p.Cl1_2 = 13.6;
    p.k01_2 = 1.3;
    p.k10_2 = p.Cl1_2 / p.V1_2;
    p.doseD2 = p.dose2 * 1e3 / p.V1_2;

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
end


function [tDose1,cDose1,tDose2,cDose2,ti,tUnion] = set_protocol(p,tmax)
    tDose1 = zeros(1,p.dosenumD1); cDose1 = zeros(1,p.dosenumD1);
    for i=1:p.dosenumD1
        tDose1(i)=p.intvlD1*(i-1)+ p.tStart1;
        cDose1(i) = p.doseD1;
    end

    tDose2 = zeros(1,p.dosenumD2); cDose2 = zeros(1,p.dosenumD2);
    for i=1:p.dosenumD2
        tDose2(i)=p.intvlD2*(i-1)+ p.tStart2;
        cDose2(i) = p.doseD2;    
    end

    % Time interval is the UNION of both tDose1 and tDose2
    tUnion = union( tDose1, tDose2, 'sorted');
    if find(tUnion==0) == 1 % do nothing, t = 0 is in there already
    else % Otherwise, include t = 0
        tUnion = [0 tUnion];
    end
    ti = [tUnion tmax];
end

function sol = solve_model(tDose1,cDose1,tDose2,cDose2,tUnion,ti,p,ICs,options) 
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

function [time, tumor] = extract_volume(sol)
        time = []; 
        x1 = []; x2 = []; x3 = []; x4 = []; 
        for j = 1:size(sol,2)
            time = [time sol{j}.x];
            x1 = [x1 sol{j}.y(11,:)]; 
            x2 = [x2 sol{j}.y(12,:)]; 
            x3 = [x3 sol{j}.y(13,:)]; 
            x4 = [x4 sol{j}.y(14,:)]; 
        end
        tumor = x1+x2+x3+x4;     
end

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

function data = load_experimental_data()
    % === Pembro PK ===
    data.pembro.time_10mpk = [ ...
        0.23 0.24 0.05 0.25 0.27 0.28 0.39 0.78 0.71 0.99 3.12 3.14 5.15 ...
        7.21 7.03 6.97 7.17 7.09 7.30 8.10 7.82 7.14 7.15 7.93 10.16 10.08 ...
        12.03 12.04 14.07 13.96 13.51 14.11 13.77 14.06 13.70 14.10 14.88 ...
        13.91 14.12 14.71 15.02 14.66 14.97 17.03 17.05 17.06 18.88 19.00 ...
        20.90 21.10 20.94 ];

    data.pembro.conc_10mpk = [ ...
        137.73 127.31 120.37 114.58 98.38 90.28 82.18 75.23 53.24 59.03 ...
        55.56 41.67 56.71 28.94 185.19 156.25 146.99 136.57 121.53 97.22 ...
        91.44 93.75 84.49 74.07 71.76 59.03 39.35 28.94 28.94 42.82 259.26 ...
        241.90 209.49 202.55 181.71 173.61 165.51 164.35 152.78 145.83 ...
        125.00 108.80 89.12 57.87 47.45 37.04 49.77 25.46 54.40 47.45 25.46 ];

    data.pembro.time_1mpk = [ ...
        0.03 0.04 0.30 0.03 0.23 0.25 0.42 0.68 0.60 0.86 1.11 3.20 3.13 ...
        5.12 5.14 5.15 7.15 7.14 7.06 7.07 7.08 7.27 7.37 7.38 7.79 7.31 ...
        7.81 8.07 8.08 8.18 10.10 10.11 10.21 12.22 12.22 14.24 14.15 ...
        13.90 14.25 14.01 13.96 14.15 14.43 14.49 14.84 15.16 15.05 15.14 ...
        17.27 16.96 19.12 19.14 21.22 21.16 ];

    data.pembro.conc_1mpk = [ ...
        15.69 14.77 14.13 13.12 11.65 10.28 8.99 8.44 7.80 6.97 6.42 5.41 ...
        4.50 4.86 3.49 1.47 2.02 3.49 23.12 21.01 19.82 16.97 15.78 14.31 ...
        14.59 12.75 12.48 11.10 10.09 8.35 7.52 6.15 4.59 3.58 2.84 0.83 ...
        1.74 21.38 19.82 19.08 14.50 12.11 7.98 11.28 8.81 10.92 14.13 ...
        12.94 7.80 4.59 5.69 2.84 3.21 0.83 ];
    

    % === Bevacizumab PK ===
    data.beva.iv_days = [
        0.913687049	0.972179311	0.971097941	1.029688509	1.058688875	...
        1.177049687	1.265820296	1.951015364	...
        1.979524197	2.277294048	2.337949049	2.992275932	2.962685729	...
        3.97524085	4.927534197	5.01699295	...
        6.02905654	9.006656742	11.9835688	14.99233211	14.99233211
    ];

    data.beva.iv_conc = [
        207.2479214	161.2488252	125.4836659	99.8836547	83.22340584	...
        69.32185114	60.4419756	61.69791428	...
        45.86879043	44.79145023	57.54667334	45.71774746	47.85514734	...
        45.57162259	37.00314449	37.84539904	...
        32.15693649	24.78290306	16.28259914	17.2651766	17.2651766
    ];

    data.beva.sc_days = [
        1.004718703	1.006094992	1.014451029	1.015532399	1.109414937	...
        1.080709491	1.143920456	1.293542679	...
        1.294820661	1.921916692	2.2796534	2.993750527	2.993160689	...
        3.947026936	3.975535769	4.989073954	...
        4.987795972	5.970957504	5.970465972	8.948950931	9.007443193	...
        11.92664944	11.98465017	14.96392158	...
        14.99282364
    ];

    data.beva.sc_conc = [
        0.305265874	0.420039538	2.916503565	3.747760876	10.69263218 ...
    	13.7415626	31.93620791	37.44369059	...
        50.36046269	72.37999878	77.41329432	64.35726656	56.12966934	...
        65.63722929	48.79744072	58.36754006	...
        43.39706178	45.27638	40.39870473	38.22545035	29.74123418	...
        30.13914807	20.92344023	23.75883114	...
        19.34974652
    ];

    % H1299 
    data.H1299.control.days   = [9 12 15 18 21 24 27 30];
    data.H1299.control.volume = [28 83 133 221 400 528 696 1070];
    data.H1299.pembro.days    = [9 12 15 18 21 24 27 30];
    data.H1299.pembro.volume  = [28 66 105 115 181 231 297 448];
    data.H1299.beva.days      = [9 12 15 18 21 24 27 30];
    data.H1299.beva.volume    = [28 66 105 143 204 343 449 684];
    data.H1299.combo.days     = [9 12 15 18 21 24 27 30];
    data.H1299.combo.volume   = [22 55 71 81 92 108 129 184]; 
end

