% function void = main_MMmodel_lambda_f_heatmap_update070522(void)

% code loglikelihood (of model given GE data )
% as function of (M,K) - create 95% confidence region

% -- updated 07/05/19 --

%%
clear all; close all; clc;

%%
tic
save_file_ans = 1;
% 0: don't save
% 1: save

filename = 'loglikelihood_M_K_070522';

fprintf('Creating loglikelihood heatmap: M vs. K ... \n\n');

%%
% Michaelis-Menten (MM) 
model_ind = 3;

% 2e6 cells
Ncells = 2e6;

% approximate number of DIPs to produce one WT
% eps = 1/25.9589;


%% mesh m, K
m_vector = linspace(500,2500,300);
K_vector = linspace(2,20,300);
% K_vector = 2:0.25:20;
% M_vector = 500:50:2500;

loglikelihood_values_MMmodel = zeros(length(m_vector),length(K_vector));


%% load MM fitting results to TCID50 data
load('./results/FittingResults_GEpersample_070122.mat');

%%
% fitting to GE data

% GE/mL (MA) after 18 hpi
data_hpi_18hrs = data.GE_output;

% WT MOI values
actual_moi = data.actual_moi;
actual_moi_array = actual_moi;

% number of data points
size_data = size(data_hpi_18hrs);
n_pts = size_data(1)*size_data(2);

% reshape data as list
data_hpi_18hrs_list = reshape(data_hpi_18hrs,[1,n_pts]);
% use coefficient of variation for noise
GE_noise = std(log10(data_hpi_18hrs'));
%%
% MLE parameter values: m,K - fitting to GE/mL vs. WT MOI
MM_pars_fit = results.nu_MM;


%%
for ii = 1:length(m_vector)
    ii
    %             length(m_vector)
    for jj = 1:length(K_vector)
        
        %         j
        ordered_pair = [m_vector(ii) K_vector(jj)];
        negloglikelihood_value = feval(@(y)ge_inputoutput_models_negloglikelihood(ordered_pair, data_hpi_18hrs, actual_moi_array, Ncells, GE_noise, model_ind));
        loglikelihood_values_MMmodel(ii,jj) = - negloglikelihood_value;
        
    end
end

%% MLEs
[max_logL_MMmodel ind_m] = max(max(loglikelihood_values_MMmodel'));
m_max_MM = m_vector(ind_m);
[max_logL_MMmodel ind_K] = max(max(loglikelihood_values_MMmodel));
K_max_MM = K_vector(ind_K);

% two degrees of freedom for chi-squared value
% for likelihood ratio test
chi_squared_val = 5.99;

max_logL = max_logL_MMmodel;
loglikelihood_values = loglikelihood_values_MMmodel;

loglikelihood_values_binary = NaN*loglikelihood_values;
critical_value = max_logL-chi_squared_val/2;
locs_positive = find(loglikelihood_values>critical_value);
locs_negative = find(loglikelihood_values<critical_value);
loglikelihood_values_binary(locs_positive) = 1;
loglikelihood_values_binary(locs_negative) = 0;


%%
hm_m8 = figure(5);
h = pcolor(K_vector,m_vector,loglikelihood_values_MMmodel); hold on;
contour(K_vector, m_vector,loglikelihood_values_binary,1,'k'); hold on;
plot(K_max_MM,m_max_MM,'k.','MarkerSize',30);
axis([K_vector(1) K_vector(end) m_vector(1) m_vector(end)])
set(h, 'EdgeColor', 'none');
xlabel('K'); ylabel('m')
colorbar
colormap('parula')
ax = gca; % current axes
ax.FontSize = 16;
ax.FontWeight = 'bold';


%% save file??
if save_file_ans
    
    % close all figures first
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
