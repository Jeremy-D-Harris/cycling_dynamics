% function void=main_simulate_passagestudy_allmodels(void)

clear all; close all; clc;

%%
save_file_ans = 0;
% 0: don't save
% 1: save

filename = 'Results_passagestudy_allmodels_070722';

fprintf('Simulating passage study (all models)... \n\n');

% light blue, orange
color_palette = [0.5216    0.7529    0.9765; 0.9608    0.4745    0.2275];

%% load fitted parameters
% load('./results/Results_TCID50_WTDIratios_070422.mat');
load('./results/Results_allmodels_TCID50_WTDIratios_070522.mat');



%% which model: 1-8

% parameters and initial conditions for all models
params.passages = 1:180;
params.num_pass = length(params.passages);
params.c = 2.9425e6;                                % 2 million host cells; given by Chris Brooke. Zwart model had c = 1e4
params.c_vals = params.c*ones(params.num_pass,1);
params.HAU_particles_factor = 50300000;             % given by Chris Brooke
params.LOD_TCID50 = 26;                             % limit of detection for TCID50 assay; given as 26 by Chris Brooke; for now keep all measurements above LOD by setting LOD to 0
params.LOD_HAU = 1;                                 % limit of detection for HAU assay; given as 1 by Chris Brooke; for now keep all measurements above LOD by setting LOD to 0
params.mu = 0;                                         % proportion of WT progeny that mutate into DIPS; Zwart model had mu = 0.78

%% models 1-2: input-independent models
nu_fit = results.nu_constant;

%% model 1: input-independent + frequency-independent
which_model=1

phi_m1 = results.pars_fit_m1(1);
f_m1 = results.pars_fit_m1(2);

% parameters and initial conditions
params.nu = nu_fit;
params.f = f_m1;
params.phi = phi_m1;
params.W_init = 3.0563e6;
params.D_init = 2.4342e8;


[W, D, TCID50, HAU] = simulate_passagestudy_models_deterministic_consolidate(params,which_model);
W_m1 = W;
D_m1 = D;
TCID50_m1 = TCID50;
HAU_m1 = HAU;

%% model 2: input-independent + frequency-dependent
which_model=2

lambda_m2 = results.pars_fit_m2(1);
f_m2 = results.pars_fit_m2(2);


% parameters and initial conditions
params.f = f_m2;
params.lambda = lambda_m2;
params.W_init = 5.4555e5;
params.D_init = 6.4088e7;

[W, D, TCID50, HAU] = simulate_passagestudy_models_deterministic_consolidate(params,which_model);
W_m2 = W;
D_m2 = D;
TCID50_m2 = TCID50;
HAU_m2 = HAU;

%% models 3-4: linear models
m_fit_linear = results.nu_linear;

%% model 3: linear + frequency-independent
which_model=3

phi_m3 = results.pars_fit_m3(1);
f_m3 = results.pars_fit_m3(2);

% parameters and initial conditions
params.m = m_fit_linear;
params.f = f_m3;
params.phi = phi_m3;
params.W_init = 1.8089e6;
params.D_init = 3.9255e7;

[W, D, TCID50, HAU] = simulate_passagestudy_models_deterministic_consolidate(params,which_model);
W_m3 = W;
D_m3 = D;
TCID50_m3 = TCID50;
HAU_m3 = HAU;

%% model 4: linear + frequency-dependent
which_model=4

lambda_m4 = results.pars_fit_m4(1);
f_m4 = results.pars_fit_m4(2);

% parameters and initial conditions
params.f = f_m4;
params.lambda = lambda_m4;
params.W_init = 2.1393e10;
params.D_init = 1.2813e11;

[W, D, TCID50, HAU] = simulate_passagestudy_models_deterministic_consolidate(params,which_model);
W_m4 = W;
D_m4 = D;
TCID50_m4 = TCID50;
HAU_m4 = HAU;

%% models 5-6: Saturating Hill models
m_fit_Hill = results.nu_Hill(1);
K_fit_Hill = results.nu_Hill(2);
n_fit = results.nu_Hill(3);

%% model 5: Hill + input-independent - need to refit params
which_model=5

phi_m5 = results.pars_fit_m5(1);
f_m5 = results.pars_fit_m5(2);


% parameters and initial conditions
params.m = m_fit_Hill;
params.K = K_fit_Hill;
params.n = n_fit;
params.f = f_m5;
params.phi = phi_m5;
params.W_init = 1.8089e6;
params.D_init = 3.9255e7;

[W, D, TCID50, HAU] = simulate_passagestudy_models_deterministic_consolidate(params,which_model);
W_m5 = W;
D_m5 = D;
TCID50_m5 = TCID50;
HAU_m5 = HAU;

%% model 6: Hill + frequency-dependent - need to refit params
which_model=6

lambda_m6 = results.pars_fit_m6(1);
f_m6 = results.pars_fit_m6(2);


% parameters and initial conditions
params.f = f_m6;
params.lambda = lambda_m6;
params.W_init = 2.1393e10;
params.D_init = 1.2813e11;

[W, D, TCID50, HAU] = simulate_passagestudy_models_deterministic_consolidate(params,which_model);
W_m6 = W;
D_m6 = D;
TCID50_m6 = TCID50;
HAU_m6 = HAU;

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
params.W_init = 3940238.98090757;
params.D_init = 247030396.745689;

[W, D, TCID50, HAU] = simulate_passagestudy_models_deterministic_consolidate(params,which_model);
W_m7 = W;
D_m7 = D;
TCID50_m7 = TCID50;
HAU_m7 = HAU;


%% model 8: MM model + frequency-dependent
which_model=8

lambda_fit = results.pars_fit_m8(1);
f_fit_m8 = results.pars_fit_m8(2);


% parameters and initial conditions
params.f = f_fit_m8;
params.lambda = lambda_fit;
params.mu = 0;
% params.W_init = 46549069.0993668;
% params.D_init = 1569845936.55439;
params.W_init = 3.0563e6;
params.D_init = 2.4342e8;

[W, D, TCID50, HAU] = simulate_passagestudy_models_deterministic_consolidate(params,which_model);
W_m8 = W;
D_m8 = D;
TCID50_m8 = TCID50;
HAU_m8 = HAU;

%% now plot models
% 4 panels: input-independent,linear,Hill,MM

figure(1);
set(gcf, 'Position',  [100, 200, 1250, 300]);

%% input-independent
figure(1); subplot(1,4,1);
p(1)=loglog(TCID50_m1, HAU_m1, '.-','Color',color_palette(1,:),'MarkerSize',15,'LineWidth',2); hold on;
p(2)=loglog(TCID50_m2, HAU_m2, '.-','Color',color_palette(2,:),'MarkerSize',15,'LineWidth',2); hold on;

xlabel('TCID50'); ylabel('Total Particles');
title('Input-independent');
axis([10^1 10^8 10^-3 10^6]);
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 16;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

%% linear
figure(1); subplot(1,4,2);
p(1)=loglog(TCID50_m3, HAU_m3, '.-','Color',color_palette(1,:),'MarkerSize',15,'LineWidth',2); hold on;
p(2)=loglog(TCID50_m4, HAU_m4, '.-','Color',color_palette(2,:),'MarkerSize',15,'LineWidth',2); hold on;

xlabel('TCID50'); ylabel('Total Particles');
title('Linear');
axis([10^1 10^8 10^-3 10^6]);
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 16;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

%% Hill
figure(1); subplot(1,4,3);
p(1)=loglog(TCID50_m5, HAU_m5, '.-','Color',color_palette(1,:),'MarkerSize',15,'LineWidth',2); hold on;
p(2)=loglog(TCID50_m6, HAU_m6, '.-','Color',color_palette(2,:),'MarkerSize',15,'LineWidth',2); hold on;

xlabel('TCID50'); ylabel('Total Particles');
title('Hill');
axis([10^1 10^8 10^-3 10^6]);
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 16;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

%% MM
figure(1); subplot(1,4,4);
p(1)=loglog(TCID50_m7, HAU_m7, '.-','Color',color_palette(1,:),'MarkerSize',15,'LineWidth',2); hold on;
p(2)=loglog(TCID50_m8, HAU_m8, '.-','Color',color_palette(2,:),'MarkerSize',15,'LineWidth',2); hold on;

xlabel('TCID50'); ylabel('Total Particles');
title('Michaelis-Menten');
axis([10^1 10^8 10^-3 10^6]);
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 16;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

%% collect results

results.HAU_m1 = HAU_m1;
results.HAU_m2 = HAU_m2;
results.HAU_m3 = HAU_m3;
results.HAU_m4 = HAU_m4;
results.HAU_m5 = HAU_m5;
results.HAU_m6 = HAU_m6;
results.HAU_m7 = HAU_m7;
results.HAU_m8 = HAU_m8;

results.TCID50_m1 = TCID50_m1;
results.TCID50_m2 = TCID50_m2;
results.TCID50_m3 = TCID50_m3;
results.TCID50_m4 = TCID50_m4;
results.TCID50_m5 = TCID50_m5;
results.TCID50_m6 = TCID50_m6;
results.TCID50_m7 = TCID50_m7;
results.TCID50_m8 = TCID50_m8;



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
