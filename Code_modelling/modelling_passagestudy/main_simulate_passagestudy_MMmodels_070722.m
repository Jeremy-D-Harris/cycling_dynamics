% function void=main_simulate_passagestudy_allmodels(void)

clear all; close all; clc;

%%
save_file_ans = 0;
% 0: don't save
% 1: save

filename = 'Results_passagestudy_MMmodels_070722';

fprintf('Simulating passage study (Michaelis-Menten forms)... \n\n');

% light blue, orange
color_palette = [0.5216    0.7529    0.9765; 0.9608    0.4745    0.2275]; 

%% load fitted parameters
% load('./results/Results_TCID50_WTDIratios_070422.mat');
load('./results/Results_MMmodels_TCID50_WTDIratios_070622.mat');

%% Michaelis-Menten models: 7-8

% parameters and initial conditions for all models
params.passages = 1:180;
params.num_pass = length(params.passages);
params.c = 2.9425e6;                                % 2 million host cells; given by Chris Brooke. Zwart model had c = 1e4
params.c_vals = params.c*ones(params.num_pass,1);
params.HAU_particles_factor = 50300000;             % given by Chris Brooke
params.LOD_TCID50 = 26;                             % limit of detection for TCID50 assay; given as 26 by Chris Brooke; for now keep all measurements above LOD by setting LOD to 0
params.LOD_HAU = 1;                                 % limit of detection for HAU assay; given as 1 by Chris Brooke; for now keep all measurements above LOD by setting LOD to 0
params.mu = 0;                                         % proportion of WT progeny that mutate into DIPS; Zwart model had mu = 0.78
% params.W_init = 1e9;
% params.D_init = 1e6;
params.W_init = 3.4204817e7;
params.D_init = 3.41042827e8;
% params.W_init = 46549069.0993668;
% params.D_init = 1569845936.55439;


%% models 7-8: Saturating MM models
m_fit_MM = results.nu_MM(1);
K_fit_MM = results.nu_MM(2);

%% model 7: MM model + frequency-independent
which_model=7

phi_fit = results.pars_fit_m7(1);
f_fit_m7 = results.pars_fit_m7(2);


% parameters and initial conditions
params.m = m_fit_MM;
params.K = K_fit_MM;
params.f = f_fit_m7;
params.phi = phi_fit;
params.mu = 0;


[W, D, TCID50, HAU] = simulate_passagestudy_models_deterministic_consolidate(params,which_model);
W_m7 = W;
D_m7 = D;
TCID50_m7 = TCID50;
HAU_m7 = HAU;


%% model 8: MM model + frequency-dependent
which_model=8

% lambda_fit = 4.5;
lambda_fit = results.pars_fit_m8(1);
f_fit_m8 = results.pars_fit_m8(2);


% parameters and initial conditions
params.f = f_fit_m8;
params.lambda = lambda_fit;
params.mu = 0;

[W, D, TCID50, HAU] = simulate_passagestudy_models_deterministic_consolidate(params,which_model);
W_m8 = W;
D_m8 = D;
TCID50_m8 = TCID50;
HAU_m8 = HAU;

%% collect results
results.W_m7 = W_m7;
results.D_m7 = D_m7;
results.TCID50_m7 = TCID50_m7;
results.HAU_m7 = HAU_m7;

results.W_m8 = W_m8;
results.D_m8 = D_m8;
results.TCID50_m8 = TCID50_m8;
results.HAU_m8 = HAU_m8;

%% now plot models 
% 3 panels: two MM models
% Particles, TCID50, Particles/TCID50

figure(1);
% set(gcf, 'Position',  [100, 200, 1250, 300]);
set(gcf, 'Position',  [100, 200, 725, 1200])

figure(1); 
subplot(3,1,1); 
p(1) = semilogy(params.passages, HAU_m7*params.HAU_particles_factor, '.-','Color',color_palette(1,:),'MarkerSize',15,'LineWidth',2);hold on;
p(2) = semilogy(params.passages, HAU_m8*params.HAU_particles_factor, '.-','Color',color_palette(2,:),'MarkerSize',15,'LineWidth',2);hold on;
axis([0 75 10^6 10^10]);
xlabel('Passage Number'); ylabel('Particles');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 16;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

subplot(3,1,2); 
semilogy(params.passages, TCID50_m7, '.-','Color',color_palette(1,:),'MarkerSize',15,'LineWidth',2);hold on
semilogy(params.passages, TCID50_m8, '.-','Color',color_palette(2,:),'MarkerSize',15,'LineWidth',2);hold on
axis([0 75 10^0 10^10]);
xlabel('Passage Number'); ylabel('TCID50');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 16;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

subplot(3,1,3); 
semilogy(params.passages, params.HAU_particles_factor*HAU_m7./TCID50_m7, '.-','Color',color_palette(1,:),'MarkerSize',15,'LineWidth',2);hold on;
semilogy(params.passages, params.HAU_particles_factor*HAU_m8./TCID50_m8, '.-','Color',color_palette(2,:),'MarkerSize',15,'LineWidth',2);hold on;
axis([0 75 10^0 10^8]);
xlabel('Passage Number'); ylabel('Particles/TCID50');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 16;
f1.FontWeight = 'normal';
f1.FontName = 'Times';


%% save file??
if save_file_ans
    
    % save for plotting
    folder_location = '../../Code_plt_ms_figures/results/';
    save(strcat(folder_location,filename),'results','params');
    
    fprintf('File saved:\n');
    fprintf(strcat(filename,'\n\n'));
    
    fprintf('Location:\n');
    fprintf(strcat(folder_location,'\n\n'));
    
else
    
    fprintf('File not saved.\n');
    
end