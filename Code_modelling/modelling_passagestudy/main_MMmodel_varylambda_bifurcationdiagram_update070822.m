% function void = main_MMmodel_varylambda_bifurcationdiagram_update070622(void)

% bifurcation diagram as function of lambda

% -- updated 07/08/19 --

%%
clear all; close all; clc;

%%
tic
save_file_ans = 0;
% 0: don't save
% 1: save

filename = 'bifurcationdiagram_lambda_070822';

fprintf('Creating bifurcation diagram: cyclying amplitude vs. lambda ... \n\n');

%%
% Michaelis-Menten (MM) + freq-dependent
which_model = 8;

num_passages = 1:280;


%% mesh lambda, f
lambda_vector = linspace(0.5,20,80);
% f_vector = linspace(0.01,0.1,200);

HAU_varylambda = zeros(length(lambda_vector),length(num_passages));
TCID50_varylambda = zeros(length(lambda_vector),length(num_passages));


%% load MM fitting results to TCID50 data
% load('./results/Results_MMmodels_TCID50_WTDIratios_070722.mat');
% load('./results/Results_MMmodels_TCID50_WTDIratios_070622.mat');
load('./results/Results_allmodels_TCID50_WTDIratios_070522.mat');


% %%
% % fitting to TCID50 data
% newdata_TCID50 = data_TCID50.newdata_TCID50;
%
% % TCID50_noise = std(log10(newdata_TCID50'));
%
% % WT MOIs (MA)
% WT_mois_MA = data_TCID50.WT_mois_MA;
%
% % DIP MOIs (MA)
% DI_mois_MA = data_TCID50.DI_mois_MA;

%%
% MLE parameter values: m,K - fitting to GE/mL vs. WT MOI
MM_pars_fit = results.nu_MM;
params.m = MM_pars_fit(1);
params.K = MM_pars_fit(2);

% MLE parameter values: lambda, f - fitting to TCID50 vs. WT:DI MOI ratios
two_pars_fit = results.pars_fit_m8;
% params.lambda = two_pars_fit(1);
params.f = two_pars_fit(2);

% parameters and initial conditions for all models
params.passages = num_passages;
params.num_pass = length(params.passages);
params.c = 2.9425e6;                                % 2 million host cells; given by Chris Brooke. Zwart model had c = 1e4
params.c_vals = params.c*ones(params.num_pass,1);
params.HAU_particles_factor = 50300000;             % given by Chris Brooke
params.LOD_TCID50 = 26;                             % limit of detection for TCID50 assay; given as 26 by Chris Brooke; for now keep all measurements above LOD by setting LOD to 0
params.LOD_HAU = 1;                                 % limit of detection for HAU assay; given as 1 by Chris Brooke; for now keep all measurements above LOD by setting LOD to 0
params.mu = 0;                                         % proportion of WT progeny that mutate into DIPS; Zwart model had mu = 0.78


%%
for ii = 1:length(lambda_vector)
    
    ii
    
    this_lambda = lambda_vector(ii);
    params.lambda = this_lambda;
    if this_lambda < 4
        
%         params.W_init = 1.811857e6;
%         params.D_init = 1.51313933e8;
        params.W_init = 3.4204817e7;
        params.D_init = 3.41042827e8;
        
    else
        params.W_init = 3.0563e6;
        params.D_init = 2.4342e8;
        
    end
    
    [this_W, this_D, this_TCID50, this_HAU] = simulate_passagestudy_models_deterministic_consolidate(params,which_model);
    HAU_varylambda(ii,:) = this_HAU;
    TCID50_varylambda(ii,:) = this_TCID50;
    max_HAU_varylambda(ii) = max(this_HAU(200:end));
    min_HAU_varylambda(ii) = min(this_HAU(200:end));
    
    
end


%% now plot the results
figure(1); set(gcf, 'Position',  [200, 100, 412.5, 350]);
semilogy(lambda_vector,max_HAU_varylambda,'k.-'); hold on;
semilogy(lambda_vector,min_HAU_varylambda,'k.-');

%% panel A
% subplot(1,2,1);

xlabel('$\lambda$','interpreter','latex'); ylabel('HAU');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';



%% save file??
if save_file_ans
    
    close all;
    
    % save for plotting
    folder_location = '../../Code_plt_ms_figures/results/';
    save(strcat(folder_location,filename));
    %     save(strcat(folder_location,filename),'-v7.3');
    
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
