%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% The user is prompted to enter what behavior they see in their         %
% "patient" in response to the following initiation protocol:           %
% 1) The initiation protocol will start by administering a standard     %
%    dose of drug 2 (beva)* and waiting tmax_d1 time to measure the     %
%    tumor volume. The user enters the relative change (compared to     %
%    baseline) for their "patient", defined as (Vf-V0)/V0 observed in   %
%    the patient "data". The user is told what an acceptable range for  %
%    this relative change is**, and if a value outside the range is     %
%    entered, the code will throw an error message and exit. Otherwise, %
%    it computes the best fit value of k2_d2 given the relative change  %
%    observed and entered for the patient.                              %
% 2) The next part of the initiation protocol is to administer a        %
%    standard dose of drug 1 (pembro), once tmax_d2 time has passed     %
%    from the adminstration of dbeva. Relative change is computed from  %
%    time this dose is given (with all ICs set based on variable        %
%    values at time tmax_d1) over a period of time tmax_d2. From this,  %
%    it computes the best fit value of k2_d1 given the relative change  %
%    observed and entered for the patient.                              %
% 3) Finally, a standard combination dose is administered once tmax_d2  %
%    time has passed from the administration of pembro. The range of    %
%    allowable responses (based on bounds of combination parameter, k3) %
%    is computed for this patient (knowing their k2_d1, k2_d2).         %
%    Relative change is computed from time this combination is given    %
%    (with all ICs set based on variable values at time tmax_d1 +       %
%    tmax_d2) over a period of time tmax_combo. From this, it computes  %
%    the best fit value of k3 given the relative change observed and    %
%    entered for the patient.                                           %
%                                                                       %
% The code also reads in data processed in analyze_data_3D_TGIthresh.m  %
% and matches the patients (based on their fit values of k2_d1, k2_d2,  %
% k3) to their near Pareto optimal doses.                               %
% - The response to the entirety of dosing space is shown in a heatmap  %
%   that includes the TGI associated with the protocol.                 %
% - The "top" combination dose is also computed, where this is defined  %
%   the combination that is most likely to be Pareto optimal across     %
%   four synergy spaces, for parameterizations that fall within the     %
%   voxel.                                                              %
% - Not only is the top combination therapy displayed, but so are all   %
%   doses for which probability of being Pareto optimal across synergy  %
%   spaces (the Pareto optimality score, for that voxel) is within      %
%   (threshold_for_opt*100)% of the optimal probability.                %
%                                                                       %   
% Notes:                                                                %
% * We administer beva first, as it has the shorter half life and thus  %
%   we have to wait less time for it to meet our clearance threshold of %
%   levels < 90% of Cmax. (This time is 24 days for beva but 31 days    %
%   for pembro.                                                         %        
% **How allowable range is calculated: Efficacy is given as change      %
%   (Vf-V0) relative to baseline (V0, reset each time a new dose is     %
%   given) in response single dose.                                     %
%                                                                       %
% Case used for paper                                                   %
% 1) Beva and wait 24 days: relative change = 70 gives                  %
%    k2_d2 = 0.129822                                                   %
% 2) Pembro and wait 6 days: relative change = 0.2 gives                %
%    k2_d1 = 0.148127                                                   %
% 3) Combination and wait 6 days: relative change = -0.3 gives          %
%    k3 = 0.040989                                                      %
%                                                                       %
% Updated 8/5/2025.                                                     %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read in data and create a directory for the patient
clear all; close all; clc;
path_read = 'Output_H1299_3D';
options = odeset('RelTol', 1.0e-10,'MaxStep',0.1); 
options_fit = optimoptions(@fmincon,'OptimalityTolerance',1e-10,...
    'Display','notify-detailed'); 

% Directory to save figures and output
cl = clock; clN = 0;
for ii = 2:5 % Month/day/hour/minute
    clN = floor(100*clN + cl(ii));
end
path = [path_read '/Patient_' , num2str(clN) '_3Doses'];
if exist(path, 'dir') ~= 7
    mkdir(path)
end
myfile = "patient_output.txt";
filepath = fullfile(path, myfile);
diary(filepath);
mat_file1 = [path_read '/output.mat'];
load(mat_file1);  

%% Let tumor grow for 9 days without drug so baseline volume is measurable
tmax = 9;
tmax_adjust = tmax;
dose1 = 0; % mg/kg; pembro  
dose2 = 0; % bevacizumab
[tDose1,cDose1,tDose2,cDose2,ti,tUnion] = set_protocol(p,tmax,dose1,dose2);

sol_no = solve_model(0,tDose1,cDose1,tDose2,cDose2,...
    tUnion,ti,p,ICs,options,0); 
time_sim_no = []; x1_no = []; x2_no = []; x3_no = []; x4_no = []; 
D1iv_no = []; D1p_no = []; D1t_no = []; D2iv_no = []; D2p_no = [];                                                            
for j = 1:size(sol_no,2)
    time_sim_no = [time_sim_no sol_no{j}.x];
    D1iv_no = [D1iv_no sol_no{j}.y(1,:)];   % D1_iv = y(1);
    D1p_no  = [D1p_no sol_no{j}.y(2,:)];    % D1_p  = y(2);
    D1t_no  = [D1t_no sol_no{j}.y(3,:)];    % D1_t  = y(3)
    D2iv_no = [D2iv_no sol_no{j}.y(6,:)];   % D2_iv = y(6);
    D2p_no  = [D2p_no sol_no{j}.y(7,:)];    % D2_p  = y(7);
    x1_no   = [x1_no sol_no{j}.y(11,:)]; 
    x2_no   = [x2_no sol_no{j}.y(12,:)]; 
    x3_no   = [x3_no sol_no{j}.y(13,:)]; 
    x4_no   = [x4_no sol_no{j}.y(14,:)]; 
end
tumor_sim_no = x1_no + x2_no + x3_no + x4_no;
tumor_sim_no_end = tumor_sim_no(end);
V0 = tumor_sim_no(end);
fprintf('Tumor has grown from volume of %f to %f mm^3 over %d days\n',...
    ICs(11),V0,tmax);

%% Set initiation protocol to match experimental dose, but only given once
tmax_d1 = 24; % time required to get drug 1 < 90% max value
tmax_d2 = 6; % time required to get drug 2 < 90% max value = 24
tmax_d1again = 6; % arbitrary wait time - no concern about clearance since
                % this is the last measurement we will take. Initially 7 
p.dose1 = 10; % mg/kg; pembro  
p.dose2 = 1; % bevacizumab
num_pts_sens = 11; % use an odd number so best fit is included in vector
sens_vary = 0.05; % percent variation from optimal parameter
threshold_for_opt = 0.25; % controls how many Pareto optimal pts print
best_params = zeros(1,size(param_vec,1));

%% Drug 2 Monotherapy: beva
ICs(11) = V0; % reflects tumor size after 9 days of growth
protocol_num = 2;
tmax = tmax_d1;
param_min = param_vec(protocol_num,1);
param_max = param_vec(protocol_num,end);
param_guess = (param_max-param_min)/2;
fprintf('Administer beva for %d days:\n',tmax);

% Must figure out min and max response achievable for combination for 
% this patient by using model to find behavior at k2_d2 extremes
dose1 = 0; dose2 = p.dose2; % Only give drug 2
[tDose1,cDose1,tDose2,cDose2,ti,tUnion] = set_protocol(p,tmax,dose1,dose2);
[time_k2d2_min,tumor_k2d2_min] = soln_only(param_min,tDose1,cDose1,...
    tDose2,cDose2,tUnion,ti,p,ICs,options,protocol_num);  
response_min = (tumor_k2d2_min(end)-V0)/V0;
fprintf('\tRelative change from baseline for this patient at min(k2d2) is: %f\n',...
    response_min);
[time_k2d2_max,tumor_k2d2_max] = soln_only(param_max,tDose1,cDose1,...
    tDose2,cDose2,tUnion,ti,p,ICs,options,protocol_num);  
response_max = (tumor_k2d2_max(end)-V0)/V0;
fprintf('\tRelative change from baseline for this patient at max(k2d2) is: %f\n',...
    response_max);

fprintf('\tChange relative to baseline will be in range [%f, %f]\n',...
    response_max,response_min);
fprintf('\t\t(Larger the change, less effective the therapy)\n\t');
prompt_d2 = "Enter change relative to baseline measured for patient: ";
response_d2 = input(prompt_d2); % Stores value answered to prompt
if (response_d2<response_max) || (response_d2>response_min)
    fprintf('\t\tResponse must be between established range for drug 2 - exiting\n');
    stop
end
V_tmax = V0 + V0*response_d2;  % convert back to volume at tmax
fprintf('\t\tModel will fit to point V(%d) = %f mm^3\n',tmax,V_tmax);

% Now find value of k2_d2 that best fits data point (tmax,V_tmax)
fun_d2 = @(z)objective(z,tmax,V_tmax,tDose1,cDose1,tDose2,cDose2,tUnion,...
    ti, p,ICs,options,protocol_num);
% fit_initial = objective(param_guess,tmax,V_tmax,tDose1,cDose1,tDose2,...
%     cDose2,tUnion,ti,p,ICs,options,protocol_num); 
[best_fit_params, best_fit, exitflag_d2, output_d2] = fmincon(fun_d2,...
    param_guess,[],[],[],[],param_min,param_max,[],options_fit);
fprintf('\tk2_d2 changed from initial value of %f to best-fit value of %f\n',...
    param_guess,best_fit_params);

% Visualize best-fit 
best_param(protocol_num) = best_fit_params; 
sol_d2 = solve_model(best_param(protocol_num),tDose1,cDose1,tDose2,cDose2,...
    tUnion,ti,p,ICs,options,protocol_num); 
time_sim_d2 = []; x1_d2 = []; x2_d2 = []; x3_d2 = []; x4_d2 = []; 
D1iv_d2 = []; D1p_d2 = []; D1t_d2 = []; D2iv_d2 = []; D2p_d2 = [];                                                            
for j = 1:size(sol_d2,2)
    time_sim_d2 = [time_sim_d2 sol_d2{j}.x];
    D1iv_d2 = [D1iv_d2 sol_d2{j}.y(1,:)];   % D1_iv = y(1);
    D1p_d2  = [D1p_d2 sol_d2{j}.y(2,:)];    % D1_p  = y(2);
    D1t_d2  = [D1t_d2 sol_d2{j}.y(3,:)];    % D1_t  = y(3)
    D2iv_d2 = [D2iv_d2 sol_d2{j}.y(6,:)];   % D2_iv = y(6);
    D2p_d2  = [D2p_d2 sol_d2{j}.y(7,:)];    % D2_p  = y(7);
    x1_d2   = [x1_d2 sol_d2{j}.y(11,:)]; 
    x2_d2   = [x2_d2 sol_d2{j}.y(12,:)]; 
    x3_d2   = [x3_d2 sol_d2{j}.y(13,:)]; 
    x4_d2   = [x4_d2 sol_d2{j}.y(14,:)]; 
end
tumor_sim_d2 = x1_d2 + x2_d2 + x3_d2 + x4_d2;
tumor_sim_d2_end = tumor_sim_d2(end);
fit_d2 = abs(tumor_sim_d2_end-V_tmax);
data_time_all = tmax;
data_tumor_all = V_tmax;

figure;
sgtitle('Beva Monotherapy','FontSize',16,'FontWeight','bold');
set(gcf, 'Units', 'Normalized','OuterPosition',[0.05, 1-0.75, 0.55, 0.75]);
subplot(2,2,1)
plot(time_sim_d2,tumor_sim_d2,'LineWidth',2); hold on;
plot(tmax,V_tmax,'o'); hold off;
xlabel('Time','FontSize',14)
ylabel('Tumor','FontSize',14)
xlim([0 tmax]);
legend('Model','Data','Location','NorthWest','FontSize',16);
title(['k2_d_2 = ' num2str(best_param(protocol_num)) ' has cost = ' ...
    num2str(fit_d2)],'FontSize',14);

% And do a quick study on how changing k2_d2 changes fit
k2_d2_vec = linspace((1-sens_vary)*best_fit_params,...
    (1+sens_vary)*best_fit_params,num_pts_sens);
fit_k2_d2 = zeros(size(k2_d2_vec));
for i = 1:length(fit_k2_d2)
    fit_k2_d2(i) = objective(k2_d2_vec(i),tmax,V_tmax,tDose1,cDose1,...
        tDose2,cDose2,tUnion,ti,p,ICs,options,protocol_num); 
end
subplot(2,2,2)
plot(k2_d2_vec,fit_k2_d2,'o-','LineWidth',2);
xlabel('k2_d_2','FontSize',14);
ylabel('Cost Function Value','FontSize',14); 
title(['Sensitivity to ' num2str(100*sens_vary) '% Variation'],...
    'FontSize',14);

subplot(2,2,3)
plot(time_sim_d2,D1p_d2,'LineWidth',2); 
xlabel('Time','FontSize',14)
ylabel('Drug 1','FontSize',14)
xlim([time_sim_d2(1) time_sim_d2(end)])
if max(D1p_d2)>0
    D1p_reduction = D1p_d2(end)/max(D1p_d1);
    title(['D1p(' num2str(time_sim_d2(end)) ')/max(D1p) = ' ...
        num2str(D1p_reduction)],'FontSize',14,'FontWeight','bold');
end

subplot(2,2,4)
plot(time_sim_d2,D2p_d2,'LineWidth',2); 
xlabel('Time','FontSize',14)
ylabel('Drug 2','FontSize',14)
xlim([time_sim_d2(1) time_sim_d2(end)])
if max(D2p_d2)>0
    D2p_reduction = D2p_d2(end)/max(D2p_d2);
    title(['D2p(' num2str(time_sim_d2(end)) ')/max(D2p) = ' ...
        num2str(D2p_reduction)],'FontSize',14,'FontWeight','bold');
end
fname_fig = [path '/k2_d2_Vf_' num2str(V_tmax)];
saveas(gcf,[fname_fig,'.fig']);
saveas(gcf,[fname_fig,'.png']);


%% Drug 1 Monotherapy: pembro
% Careful to preserve info from prior administration
ICs(1) = D1iv_d2(end); % D1_iv  = y(1);
ICs(2) = D1p_d2(end);  % D1_p  = y(2);
ICs(3) = D1t_d2(end);  % D1_t  = y(3);
ICs(6) = D2iv_d2(end); % D2_iv = y(6);
ICs(7) = D2p_d2(end);  % D2_p  = y(7);

ICs(11) = x1_d2(end);
ICs(12) = x2_d2(end);
ICs(13) = x3_d2(end);
ICs(14) = x4_d2(end);
V0 = tumor_sim_d2(end); % recompute based on current tumor size
p.k2_d2 = best_param(protocol_num); % now we know k2_d2, so fix it

protocol_num = 1;
tmax = tmax_d2;
param_min = param_vec(protocol_num,1);
param_max = param_vec(protocol_num,end);
param_guess = (param_max-param_min)/2;
fprintf('\nAdminister pembro for %d days:\n',tmax);

% Must figure out min and max response achievable for combination for 
% this patient by using model to find behavior at k2_d1 extremes
dose1 = p.dose1; dose2 = 0; % Only give drug 1
[tDose1,cDose1,tDose2,cDose2,ti,tUnion] = set_protocol(p,tmax,dose1,dose2);
[time_k2d1_min,tumor_k2d1_min] = soln_only(param_min,tDose1,cDose1,...
    tDose2,cDose2,tUnion,ti,p,ICs,options,protocol_num); 
response_min = (tumor_k2d1_min(end)-V0)/V0;
fprintf('\tRelative change from baseline for this patient at min(k2d1) is: %f\n',...
    response_min);
[time_k2d1_max,tumor_k2d1_max] = soln_only(param_max,tDose1,cDose1,...
    tDose2,cDose2,tUnion,ti,p,ICs,options,protocol_num); 
response_max = (tumor_k2d1_max(end)-V0)/V0;
fprintf('\tRelative change from baseline for this patient at max(k2d1) is: %f\n',...
    response_max);

fprintf('\tChange relative to baseline will be in range [%f, %f]\n',...
    response_max,response_min);
fprintf('\t\t(Larger the change, less effective the therapy)\n\t');
prompt_d1 = "Enter change relative to baseline measured for patient: ";
response_d1 = input(prompt_d1); % Stores value answered to prompt
if (response_d1<response_max) || (response_d1>response_min)
    fprintf('\t\tResponse must be between established range for drug 1 - exiting\n');
    stop
end
V_tmax = V0 + V0*response_d1;  % convert back to volume at tmax
fprintf('\t\tModel will fit to point V(%d) = %f mm^3\n',tmax,V_tmax);
% Now find value of k2_d1 that best fits data point (tmax,V_tmax)
fun_d1 = @(z)objective(z,tmax,V_tmax,tDose1,cDose1,tDose2,cDose2,tUnion,ti,...
    p,ICs,options,protocol_num);
[best_fit_params, best_fit, exitflag_d1, output_d1] = fmincon(fun_d1,param_guess,...
    [],[],[],[],param_min,param_max,[],options_fit);
fprintf('\tk2_d1 changed from initial value of %f to best-fit value of %f\n',...
    param_guess,best_fit_params);

% Visualize best-fit 
best_param(protocol_num) = best_fit_params; 
sol_d1 = solve_model(best_param(protocol_num),tDose1,cDose1,tDose2,cDose2,...
    tUnion,ti,p,ICs,options,protocol_num); 
time_sim_d1 = []; x1_d1 = []; x2_d1 = []; x3_d1 = []; x4_d1 = []; 
D1iv_d1 = []; D1p_d1 = []; D1t_d1 = []; D2iv_d1 = []; D2p_d1 = [];                                                            
for j = 1:size(sol_d1,2)
    time_sim_d1 = [time_sim_d1 sol_d1{j}.x];
    D1iv_d1 = [D1iv_d1 sol_d1{j}.y(1,:)];   % D1_iv = y(1);
    D1p_d1  = [D1p_d1 sol_d1{j}.y(2,:)];    % D1_p  = y(2);
    D1t_d1  = [D1t_d1 sol_d1{j}.y(3,:)];    % D1_t  = y(3)
    D2iv_d1 = [D2iv_d1 sol_d1{j}.y(6,:)];   % D2_iv = y(6);
    D2p_d1  = [D2p_d1 sol_d1{j}.y(7,:)];    % D2_p  = y(7);
    x1_d1   = [x1_d1 sol_d1{j}.y(11,:)]; 
    x2_d1   = [x2_d1 sol_d1{j}.y(12,:)]; 
    x3_d1   = [x3_d1 sol_d1{j}.y(13,:)]; 
    x4_d1   = [x4_d1 sol_d1{j}.y(14,:)]; 
end
tumor_sim_d1 = x1_d1 + x2_d1 + x3_d1 + x4_d1;
tumor_sim_d1_end = tumor_sim_d1(end);
fit_d1 = abs(tumor_sim_d1_end-V_tmax);
data_time_all = [data_time_all tmax+time_sim_d2(end)];
data_tumor_all = [data_tumor_all V_tmax];


figure;
sgtitle('Pembro Monotherapy','FontSize',16,'FontWeight','bold');
set(gcf, 'Units', 'Normalized','OuterPosition',[0.05, 1-0.75, 0.55, 0.75]);
subplot(2,2,1)
plot(time_sim_d1,tumor_sim_d1,'LineWidth',2);hold on;
plot(tmax,V_tmax,'o');  hold off;
xlabel('Time','FontSize',14)
ylabel('Tumor','FontSize',14)
xlim([0 tmax]);
legend('Model','Data','Location','NorthWest','FontSize',16);
title(['k2_d_1 = ' num2str(best_param(protocol_num)) ' has cost = ' ...
    num2str(fit_d1)],'FontSize',14);

% And do a quick study on how changing k2_d1 changes fit
k2_d1_vec = linspace((1-sens_vary)*best_fit_params,...
    (1+sens_vary)*best_fit_params,num_pts_sens);
fit_k2_d1 = zeros(size(k2_d1_vec));
for i = 1:length(fit_k2_d1)
    fit_k2_d1(i) = objective(k2_d1_vec(i),tmax,V_tmax,tDose1,cDose1,...
        tDose2,cDose2,tUnion,ti,p,ICs,options,protocol_num); 
end
subplot(2,2,2)
plot(k2_d1_vec,fit_k2_d1,'o-','LineWidth',2);
xlabel('k2_d_1','FontSize',14);
ylabel('Cost Function Value','FontSize',14); 
title(['Sensitivity to ' num2str(100*sens_vary) '% Variation'],...
    'FontSize',14);

subplot(2,2,3)
plot(time_sim_d1,D1p_d1,'LineWidth',2); 
xlabel('Time','FontSize',14)
ylabel('Drug 1','FontSize',14)
xlim([time_sim_d1(1) time_sim_d1(end)])
if max(D1p_d1)>0
    D1p_reduction = D1p_d1(end)/max(D1p_d1);
    title(['D1p(' num2str(time_sim_d1(end)) ')/max(D1p) = ' ...
        num2str(D1p_reduction)],'FontSize',14,'FontWeight','bold');
end

subplot(2,2,4)
plot(time_sim_d1,D2p_d1,'LineWidth',2); 
xlabel('Time','FontSize',14)
ylabel('Drug 2','FontSize',14)
xlim([time_sim_d1(1) time_sim_d1(end)])
if max(D2p_d1)>0
    D2p_reduction = D2p_d1(end)/max(D2p_d1);
    title(['D2p(' num2str(time_sim_d1(end)) ')/max(D2p) = ' ...
        num2str(D2p_reduction)],'FontSize',14,'FontWeight','bold');
end

fname_fig = [path '/k2_d1_Vf_' num2str(V_tmax)];
saveas(gcf,[fname_fig,'.fig']);
saveas(gcf,[fname_fig,'.png']);


%% Combination
% Careful to preserve info from prior administration
ICs(1) = D1iv_d1(end); % D1_iv  = y(1);
ICs(2) = D1p_d1(end);  % D1_p  = y(2);
ICs(3) = D1t_d1(end);  % D1_t  = y(3);
ICs(6) = D2iv_d1(end); % D2_iv = y(6);
ICs(7) = D2p_d1(end);  % D2_p  = y(7);

ICs(11) = x1_d1(end);
ICs(12) = x2_d1(end);
ICs(13) = x3_d1(end);
ICs(14) = x4_d1(end);
V0 = tumor_sim_d1(end); % recompute based on current tumor size
p.k2_d1 = best_param(protocol_num); % now we know k2_d2, so fix it

protocol_num = 3;
tmax = tmax_d1again;
param_min = param_vec(protocol_num,1);
param_max = param_vec(protocol_num,end);
param_guess = (param_max-param_min)/2;
fprintf('\nAdminister combination for %d days:\n',tmax);

% Must figure out min and max response achievable for combination for 
% this patient by using model to find behavior at k3 extremes
dose1 = p.dose1; dose2 = p.dose2; % Give both drugs
[tDose1,cDose1,tDose2,cDose2,ti,tUnion] = set_protocol(p,tmax,dose1,dose2);
[time_k3_min,tumor_k3_min] = soln_only(param_min,tDose1,cDose1,...
    tDose2,cDose2,tUnion,ti,p,ICs,options,protocol_num);   
response_min = (tumor_k3_min(end)-V0)/V0;
fprintf('\tRelative change from baseline for this patient at min(k3) is: %f\n',...
    response_min);
[time_k3_max,tumor_k3_max] = soln_only(param_max,tDose1,cDose1,...
    tDose2,cDose2,tUnion,ti,p,ICs,options,protocol_num); 
response_max = (tumor_k3_max(end)-V0)/V0;
fprintf('\tRelative change from baseline for this patient at max(k3) is: %f\n',...
    response_max);

fprintf('\tChange relative to baseline will be in range [%f, %f]\n',...
    response_max,response_min);
fprintf('\t\t(Larger the change, less effective the therapy)\n\t');
prompt_combo = "Enter change relative to baseline measured for patient: ";
response_combo = input(prompt_combo); % Stores value answered to prompt
if (response_combo<response_max) || (response_combo>response_min)
    fprintf('Response must be between established range for combo - exiting\n');
    stop
end
V_tmax = V0 + V0*response_combo;  % convert back to volume at tmax
fprintf('\t\tModel will fit to point V(%d) = %f mm^3\n',tmax,V_tmax);

% Now find value of k3 that best fits data point (tmax,V_tmax)
fun_combo = @(z)objective(z,tmax,V_tmax,tDose1,cDose1,tDose2,cDose2,tUnion,ti,...
    p,ICs,options,protocol_num);
[best_fit_params, best_fit, exitflag_combo, output_combo] = fmincon(fun_combo,param_guess,...
    [],[],[],[],param_min,param_max,[],options_fit);
fprintf('\tk3 changed from initial value of %f to best-fit value of %f\n',...
    param_guess,best_fit_params);

% Visualize best-fit 
best_param(protocol_num) = best_fit_params; 
sol_combo = solve_model(best_param(protocol_num),tDose1,cDose1,tDose2,cDose2,...
    tUnion,ti,p,ICs,options,protocol_num); 
time_sim_combo = []; x1_combo = []; x2_combo = []; x3_combo = []; x4_combo = []; 
D1iv_combo = []; D1p_combo = []; D1t_combo = []; D2iv_combo = []; D2p_combo = [];                                                            
for j = 1:size(sol_combo,2)
    time_sim_combo = [time_sim_combo sol_combo{j}.x];
    D1iv_combo = [D1iv_combo sol_combo{j}.y(1,:)];   % D1_iv = y(1);
    D1p_combo  = [D1p_combo sol_combo{j}.y(2,:)];    % D1_p  = y(2);
    D1t_combo  = [D1t_combo sol_combo{j}.y(3,:)];    % D1_t  = y(3)
    D2iv_combo = [D2iv_combo sol_combo{j}.y(6,:)];   % D2_iv = y(6);
    D2p_combo  = [D2p_combo sol_combo{j}.y(7,:)];    % D2_p  = y(7);
    x1_combo   = [x1_combo sol_combo{j}.y(11,:)]; 
    x2_combo   = [x2_combo sol_combo{j}.y(12,:)]; 
    x3_combo   = [x3_combo sol_combo{j}.y(13,:)]; 
    x4_combo   = [x4_combo sol_combo{j}.y(14,:)]; 
end
tumor_sim_combo = x1_combo + x2_combo + x3_combo + x4_combo;
tumor_sim_combo_end = tumor_sim_combo(end);
fit_combo = abs(tumor_sim_combo_end-V_tmax);
data_time_all = [data_time_all tmax+time_sim_d2(end)+time_sim_d1(end)];
%data_time_all = data_time_all + tmax_adjust; % adjust for no drug growth period
data_tumor_all = [data_tumor_all V_tmax];

figure;
sgtitle('Combination Therapy','FontSize',16,'FontWeight','bold');
set(gcf, 'Units', 'Normalized','OuterPosition',[0.05, 1-0.75, 0.55, 0.75]);
subplot(2,2,1)
plot(time_sim_combo,tumor_sim_combo,'LineWidth',2); hold on;
plot(tmax,V_tmax,'o'); hold off;
xlabel('Time','FontSize',14)
ylabel('Tumor','FontSize',14)
xlim([0 tmax]);
legend('Model','Data','Location','NorthWest','FontSize',16);
title(['k_3 = ' num2str(best_param(protocol_num)) ' has cost = ' ...
    num2str(fit_combo)],'FontSize',14);

%% And do a quick study on how changing k3 changes fit
k3_vec = linspace((1-sens_vary)*best_fit_params,...
    (1+sens_vary)*best_fit_params,num_pts_sens);
fit_k3 = zeros(size(k3_vec));
for i = 1:length(fit_k3)
    fit_k3(i) = objective(k3_vec(i),tmax,V_tmax,tDose1,cDose1,...
        tDose2,cDose2,tUnion,ti,p,ICs,options,protocol_num); 
end
subplot(2,2,2)
plot(k3_vec,fit_k3,'o-','LineWidth',2);
xlabel('k_3','FontSize',14);
ylabel('Cost Function Value','FontSize',14); 
title(['Sensitivity to ' num2str(100*sens_vary) '% Variation'],...
    'FontSize',14);

subplot(2,2,3)
plot(time_sim_combo,D1p_combo,'LineWidth',2); 
xlabel('Time','FontSize',14)
ylabel('Drug 1','FontSize',14)
xlim([time_sim_combo(1) time_sim_combo(end)])
if max(D1p_combo)>0
    D1p_reduction = D1p_combo(end)/max(D1p_combo);
    title(['D1p(' num2str(time_sim_combo(end)) ')/max(D1p) = ' ...
        num2str(D1p_reduction)],'FontSize',14,'FontWeight','bold');
end

subplot(2,2,4)
plot(time_sim_combo,D2p_combo,'LineWidth',2); 
xlabel('Time','FontSize',14)
ylabel('Drug 2','FontSize',14)
xlim([time_sim_combo(1) time_sim_combo(end)])
if max(D2p_combo)>0
    D2p_reduction = D2p_combo(end)/max(D2p_combo);
    title(['D2p(' num2str(time_sim_combo(end)) ')/max(D2p) = ' ...
        num2str(D2p_reduction)],'FontSize',14,'FontWeight','bold');
end

fname_fig = [path '/k3_Vf_' num2str(V_tmax)];
saveas(gcf,[fname_fig,'.fig']);
saveas(gcf,[fname_fig,'.png']);


%% Create summary plot
time_all = [time_sim_d2(2:end) ...
    time_sim_d1(2:end)+time_sim_d2(end) ... 
    time_sim_combo(2:end)+time_sim_d2(end)+time_sim_d1(end)];
D1p_all = [D1p_d2(2:end) D1p_d1(2:end) D1p_combo(2:end)];
D2p_all = [D2p_d2(2:end) D2p_d1(2:end) D2p_combo(2:end)];
tumor_all = [tumor_sim_d2(2:end) tumor_sim_d1(2:end) ...
    tumor_sim_combo(2:end)];

figure; 
sgtitle('Initiation Protocol','FontSize',16,'FontWeight','bold');
set(gcf, 'Units', 'Normalized','OuterPosition',[0.05, 1-0.8, 0.55, 0.8]);
subplot(2,3,1) % Drug 1 - concatinate time and profile
plot(time_all,D1p_all,'LineWidth',2); 
xlim([time_all(1) time_all(end)])
xlabel('Time (Days)','FontSize',14)
ylabel('Pembrolizumab concentration','FontSize',14)

subplot(2,3,2) % Drug 2
plot(time_all,D2p_all,'LineWidth',2); 
xlim([time_all(1) time_all(end)])
xlabel('Time (Days)','FontSize',14)
ylabel('Bevacizumab concentration','FontSize',14)

subplot(2,3,3) % Tumor
plot(time_all,tumor_all,'LineWidth',2); 
hold on;
plot(data_time_all,data_tumor_all,'o','LineWidth',2);  hold off;
xlim([time_all(1) time_all(end)])
xlabel('Time (Days)','FontSize',14)
ylabel('Tumor Volume (mm^3)','FontSize',14)
legend('Model','Data','Location','NorthWest','FontSize',16);

subplot(2,3,4) % k2_d1 profile
plot(k2_d1_vec,fit_k2_d1,'o-','LineWidth',2);
xlim([k2_d1_vec(1) k2_d1_vec(end)])
xlabel('k2_d_1','FontSize',14);
ylabel('Cost Function Value','FontSize',14); 

subplot(2,3,5) % k2_d2 profile
plot(k2_d2_vec,fit_k2_d2,'o-','LineWidth',2);
xlim([k2_d2_vec(1) k2_d2_vec(end)])
xlabel('k2_d_2','FontSize',14);
ylabel('Cost Function Value','FontSize',14); 

subplot(2,3,6) % k3_profile
plot(k3_vec,fit_k3,'o-','LineWidth',2);
xlim([k3_vec(1) k3_vec(end)])
xlabel('k_3','FontSize',14);
ylabel('Cost Function Value','FontSize',14); 
fname_fig = [path '/initiation_protocol'];
saveas(gcf,[fname_fig,'.fig']);
saveas(gcf,[fname_fig,'.png']);


%% Now find the parameter voxel associated with this patient 
format shortG
mat_file2 = [path_read '/output_processed_data.mat'];
load(mat_file2);  
num_pts = length(d1);
step_d1 = 0.5*(max(d1)-min(d1))/num_pts; 
step_d2 = 0.5*(max(d2)-min(d2))/num_pts; 

fprintf('\nClassify the patient by its three drug-response parameters:\n');
param1 = find(p1(1:end-1) <= best_param(1) & best_param(1) < p1(2:end), 1);
if ~isempty(param1)
    fprintf('\tk2_d1 = %g is between p1(%d) = %f and p1(%d) = %f\n', ...
        best_param(1), param1, p1(param1), param1+1, p1(param1+1));
else
    disp('k2_d1 = %g is outside the range of p1 - exiting\n',best_param(1));
    disp(p1);
    stop
end

param2 = find(p2(1:end-1) <= best_param(2) & best_param(2) < p2(2:end), 1);
if ~isempty(param2)
    fprintf('\tk2_d2 = %g is between p2(%d) = %f and p2(%d) = %f\n', ...
        best_param(2), param2, p2(param2), param2+1, p2(param2+1));
else
    disp('k2_d2 = %g is outside the range of p2 - exiting\n',best_param(2));
    disp(p2);
    stop
end

param3 = find(p3(1:end-1) <= best_param(3) & best_param(3) < p3(2:end), 1);
if ~isempty(param3)
    fprintf('\tk3 = %g is between p3(%d) = %f and p3(%d) = %f\n', ...
        best_param(3), param3, p3(param3), param3+1, p3(param3+1));
else
    disp('k3 = %g is outside the range of p3 - exiting\n',best_param(3));
    disp(p3);
    stop
end

% Convert patient parameterization into correct location in dose_matrices_all_pixel
param_location = zeros(length(p1_avg_vec),length(p2_avg_vec),length(p3_avg_vec));
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

fprintf('\nTop combination doses for patient with params in above range:\n');
fprintf('(Defined as all doses for which probability of being Pareto optimal\n');
fprintf('across synergy spaces is within %f%% of optimal probability)\n',...
    100*threshold_for_opt)
for i = 1:length(num_opt_sort)
    fprintf('\t%d. Drug 1 dose = %f mg/kg, drug 2 dose = %f mg/kg:\n',...
        i,d1_opt_sort(i),d2_opt_sort(i));
    fprintf('\t\tProb(Pareto optimal) = %f\n',num_opt_sort(i)/(Ncriterion*Nsamples));
    x = d1_opt_idx_sort(i);
    y = d2_opt_idx_sort(i);
    fprintf('\t\tTGI(combo) = %f\n',TGI_pixel_avg{param_pos}(x,y));
end

% And visualize all
figure; hold on;
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.35, 0.65]);
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
[M,c] = contour(d1,d2,TGI_pixel_avg{param_pos}',[0.9 0.9],'LineColor','r'); 
c.LineWidth = 3;
formatSpec = '%3.2f';
font = 'Arial';
for m = 1:length(d1)
    for n = 1:length(d2)
       r = num2str(near_dose_matrices_all_pixel{param_pos}(m,n)/(Ncriterion*Nsamples),formatSpec); 
       text(d1(m)-0.25*step_d1,d2(n)-0.1*step_d2,r,'fontname',font,...
           'color','black','fontweight','bold','FontSize',11);
    end
end
hold off; 
fname_fig = [path '/pareto_k2d1_' num2str(param1) '_k2d2_' num2str(param2)...
    '_k3_' num2str(param3)];
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png'])

diary off;
fsave = [path '/pareto_for_patient.mat']; 
save(fsave,'num_pts_sens','sens_vary','threshold_for_opt',...
    'response_d1','response_d2','response_combo','best_param',...
    'k2_d1_vec','fit_k2_d1','k2_d2_vec','fit_k2_d2','tumor_k3_min',...
    'tumor_k3_max','k3_vec','fit_k3','param1','param2','param3',...
    'param_pos','near_dose_matrices_all_pixel','TGI_pixel_avg');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Functions %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    p.dosenumD1 = 1;   % number of doses
    p.freq1=1; % day
    p.intvlD1 = d1_spacing*p.freq1; % dosing interval, days
    p.tStart1 = 0; % drug 1, time of start - set to time of first measurement in data
    d1_on_off = ones(1,p.dosenumD1); % 1 = all drug on by default

    % Drug 2 experimental protocol - bevacizumab, 1mpk
    p.dosenumD2 = 1;   % number of doses
    p.freq2=1; % day
    p.intvlD2 = d2_spacing*p.freq2; % dosing interval, days
    p.tStart2 = 0; % drug 2, time of start
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
    ti,p,ICs,options,drug_case) 

    if drug_case == 1
        p.k2_d1  = VP_params(1);
        %fprintf('\tUsing k2_d1 = %f\n',p.k2_d1);
    elseif drug_case == 2
        p.k2_d2  = VP_params(1);
        %fprintf('\tUsing k2_d2 = %f\n',p.k2_d2);
    elseif drug_case == 3
        p.k3 = VP_params(1);
        %fprintf('\tUsing k3 = %f\n',p.k3);
    end

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
        %tSpan = [ti(i) ti(i+1)];
        tSpan = ti(1):0.1:ti(i+1);
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

%% Fit function (to be minimized)
function fit = objective(params_to_fit,time_data,tumor_data,...
    tDose1,cDose1,tDose2,cDose2,tUnion,ti,p,ICs,options,drug_case)

    sol = solve_model(params_to_fit,tDose1,cDose1,tDose2,cDose2,...
        tUnion,ti,p,ICs,options,drug_case); 

    %% Extract relevant model output for fitting at "experimental" time points
    time_sim = []; x1 = []; x2 = []; x3 = []; x4 = []; 
    for j = 1:size(sol,2)
        time_sim = [time_sim sol{j}.x];
        x1 = [x1 sol{j}.y(11,:)]; 
        x2 = [x2 sol{j}.y(12,:)]; 
        x3 = [x3 sol{j}.y(13,:)]; 
        x4 = [x4 sol{j}.y(14,:)]; 
    end
    tumor_sim = x1+x2+x3+x4;
    tumor_sim_end = tumor_sim(end);
    fit = abs(tumor_sim_end-tumor_data);

    %fprintf('k2_d1 = %f gives SSE = %f\n',params_to_fit,SSE)
end

%% Model solution
function [time_sim,tumor_sim] = soln_only(params,tDose1,cDose1,...
    tDose2,cDose2,tUnion,ti,p,ICs,options,drug_case)

    sol = solve_model(params,tDose1,cDose1,tDose2,cDose2,...
        tUnion,ti,p,ICs,options,drug_case); 

    %% Extract relevant model output for fitting at "experimental" time points
    time_sim = []; x1 = []; x2 = []; x3 = []; x4 = []; 
    for j = 1:size(sol,2)
        time_sim = [time_sim sol{j}.x];
        x1 = [x1 sol{j}.y(11,:)]; 
        x2 = [x2 sol{j}.y(12,:)]; 
        x3 = [x3 sol{j}.y(13,:)]; 
        x4 = [x4 sol{j}.y(14,:)]; 
    end
    tumor_sim = x1+x2+x3+x4;
end

