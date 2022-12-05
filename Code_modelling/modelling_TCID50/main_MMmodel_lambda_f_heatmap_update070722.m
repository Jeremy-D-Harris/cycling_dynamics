% function void = main_MMmodel_lambda_f_heatmap_update070522(void)

% code loglikelihood (of model given TCID50 data using PB1 MOI)
% as function of (lambda,f) - create 95% confidence region

% -- updated 07/06/19 --

%%
clear all; close all; clc;

%%
tic
save_file_ans = 1;
% 0: don't save
% 1: save

filename = 'loglikelihood_lambda_f_070722';

fprintf('Creating loglikelihood heatmap: lambda vs. f ... \n\n');

%%
% Michaelis-Menten (MM) + freq-dependent
model_ind = 8;

% 2e6 cells
Ncells = 2e6;

% approximate number of DIPs to produce one WT
eps = 1/25.9589;


%% mesh lambda, f
lambda_vector = linspace(0.5,25,300);
f_vector = linspace(0.01,0.1,300);

loglikelihood_values_model8 = zeros(length(lambda_vector),length(f_vector));


%% load MM fitting results to TCID50 data
load('./results/Results_MMmodels_TCID50_WTDIratios_070522.mat');

%%
% fitting to TCID50 data
newdata_TCID50 = data_TCID50.newdata_TCID50;

TCID50_noise = std(log10(newdata_TCID50'));

% WT MOIs (MA)
WT_mois_MA = data_TCID50.WT_mois_MA;

% DIP MOIs (MA)
DI_mois_MA = data_TCID50.DI_mois_MA;

%%
% MLE parameter values: m,K - fitting to GE/mL vs. WT MOI
MM_pars_fit = results.nu_MM;

% MLE parameter values: lambda, f - fitting to TCID50 vs. WT:DI MOI ratios
two_pars_fit = results.pars_fit_m8;

%%
for ii = 1:length(lambda_vector)
    ii
    %             length(lambda_vector)
    for jj = 1:length(f_vector)
        
        %         j
        ordered_pair = [lambda_vector(ii) f_vector(jj)];
        negloglikelihood_value = feval(@(y)TCID50_inputoutput_models78_MOI_negloglikelihood_wDIP(ordered_pair, newdata_TCID50, TCID50_noise, WT_mois_MA, DI_mois_MA, MM_pars_fit, eps,Ncells,model_ind));
        loglikelihood_values_model8(ii,jj) = - negloglikelihood_value;
        
    end
end

%% MLEs
[max_logL_model8 ind_lambda] = max(max(loglikelihood_values_model8'));
lambda_max_m8 = lambda_vector(ind_lambda);
[max_logL_model8 ind_f_model8] = max(max(loglikelihood_values_model8));
f_max_m8 = f_vector(ind_f_model8);

% two degrees of freedom for chi-squared value
% for likelihood ratio test
chi_squared_val = 5.99;

max_logL = max_logL_model8;
loglikelihood_values = loglikelihood_values_model8;

loglikelihood_values_binary = NaN*loglikelihood_values;
critical_value = max_logL-chi_squared_val/2;
locs_positive = find(loglikelihood_values>critical_value);
locs_negative = find(loglikelihood_values<critical_value);
loglikelihood_values_binary(locs_positive) = 1;
loglikelihood_values_binary(locs_negative) = 0;


%%
hm_m8 = figure(5);
h = pcolor(f_vector,lambda_vector,loglikelihood_values_model8); hold on;
contour(f_vector, lambda_vector,loglikelihood_values_binary,1,'k'); hold on;
plot(f_max_m8,lambda_max_m8,'k.','MarkerSize',30);
axis([f_vector(1) f_vector(end) lambda_vector(1) lambda_vector(end)])
set(h, 'EdgeColor', 'none');
xlabel('f'); ylabel('lambda');
colorbar
colormap('parula')
ax = gca; % current axes
ax.FontSize = 16;
ax.FontWeight = 'bold';



%% save file??
if save_file_ans
    
    close all;
    
    % save for plotting
    folder_location = '../../Code_plt_ms_figures/results/';
    %     save(strcat(folder_location,filename),'results','data_TCID50');
    save(strcat(folder_location,filename),'-v7.3');
    
    fprintf('File saved:\n');
    fprintf(strcat(filename,'\n\n'));
    
    fprintf('Location:\n');
    fprintf(strcat(folder_location,'\n\n'));
    
    %     % save for further analysis
    %     folder_location = './results/';
    %     save(strcat(folder_location,filename),'-v7.3');
    % %     save(strcat(folder_location,filename),'results','data_TCID50');
    %
    %     fprintf('File saved:\n');
    %     fprintf(strcat(filename,'\n\n'));
    %
    %     fprintf('Location:\n');
    %     fprintf(strcat(folder_location,'\n\n'));
    
else
    
    fprintf('File not saved.\n');
    
end

toc
