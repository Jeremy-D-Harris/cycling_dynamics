% function void = main_estimate_GEpersample_050222(void)

%% estimate viral yield as a function of WT MOI
% -- updated 07/01/22 --

clear all; close all; clc;

%% save figure?
save_file_ans = 0;
% 0: don't save
% 1: save

filename = 'FittingResults_GEpersample_070122';

fprintf('Fitting models to data... \n\n');

% assume 2 million cells
Ncells = 2e6;

% violoet, orange, light blue
colors_rgb_values = [0.6627 0.3529 0.6314; 0.9608  0.4745   0.2275; 0.5216    0.7529    0.9765];
% colors_rgb_values = [0.4940, 0.1840, 0.5560; 0.8500, 0.3250, 0.0980; 0.4660, 0.6740, 0.1880];

% load data - GE/mL vs. WT MOI (MA based)
load('./data/PR8_WTonly_update091219_MAbased.mat'); 


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

% for plotting purposes - GE vs. WT MOI
actual_moi_finer = transpose(0:.01:10); 
results.actual_moi_finer =actual_moi_finer;


%% Input-independent model
which_model = 1;

fprintf('Input-independent model: \n');

x0_constant = 500; % nu_init
k_constant = length(x0_constant);
% options = optimset('MaxFunEvals',2e5,'MaxIter',2e5,'TolFun',1e-8, 'TolX',1e-8,'Display','iter');
[x, fval_constant, exitflag, output] = fminsearch(@(x)ge_inputoutput_models_negloglikelihood(x, data_hpi_18hrs, actual_moi_array, Ncells, GE_noise, which_model), x0_constant);
nu_MLE_constant = x;
logL_constant = -fval_constant;
AIC_constant = 2*(k_constant-logL_constant);

results.nu_constant=nu_MLE_constant;
results.logL_constant=logL_constant;
results.AIC_constant=AIC_constant;

fprintf('nu_MLE = %4.2f \n',nu_MLE_constant);
fprintf('Loglikelihood = %4.2f \n',logL_constant);
fprintf('AIC = %4.2f \n',AIC_constant);

% get GE constant model - best fit nu
GE_model_constant = Get_GE_models(actual_moi_finer, nu_MLE_constant, Ncells, which_model);

% rgb_values = [0.4940, 0.1840, 0.5560];
ilist = 0:15; 
cell_output_MLE_constant = ones(size(ilist))*nu_MLE_constant; 
cell_output_MLE_constant(1) = 0;

results.GE_model_constant=GE_model_constant;
results.cell_output_MLE_constant=cell_output_MLE_constant;


%% linear model:
which_model = 2;

fprintf('-------------------------------------- \n\n');
fprintf('Linear model \n\n');
      
x0_linear = 80; % slope m_init
k_linear = length(x0_linear);
% options = optimset('MaxFunEvals',2e3,'MaxIter',2e3);
% options=[];
[x, fval_linear, exitflag, output] = fminsearch(@(x)ge_inputoutput_models_negloglikelihood(x, data_hpi_18hrs, actual_moi_array, Ncells, GE_noise, which_model), x0_linear);
nu_MLE_linear = x;
logL_linear = -fval_linear;
AIC_linear = 2*(k_linear-logL_linear);

results.nu_linear=nu_MLE_linear;
results.logL_linear=logL_linear;
results.AIC_linear=AIC_linear;

fprintf('m_MLE = %4.2f \n',nu_MLE_linear);
fprintf('Loglikelihood = %4.2f \n',logL_linear);
fprintf('AIC = %4.2f \n',AIC_linear);

% get GE linear model - best fit slope, m
GE_model_linear = Get_GE_models(actual_moi_finer, nu_MLE_linear, Ncells, which_model);
cell_output_MLE_linear = nu_MLE_linear*ilist; 
% cell_output_MLE_linear(1) = 0;

results.GE_model_linear=GE_model_linear;
results.cell_output_MLE_linear=cell_output_MLE_linear;


%% Michaelis-Menten model:
which_model = 3;

fprintf('-------------------------------------- \n\n');
fprintf('Michaelis-Menten (Saturating) model \n\n');
       
x0_MM = [783.3832, 4.4587]; % MM function parameters m_init and K_init
k_MM = length(x0_MM);
% options = optimset('MaxFunEvals',2e3,'MaxIter',2e3);
% options=[];
[x, fval, exitflag, output] = fminsearch(@(x)ge_inputoutput_models_negloglikelihood(x, data_hpi_18hrs, actual_moi_array, Ncells, GE_noise, which_model), x0_MM);
m_MLE_MM = x(1);
K_MLE_MM = x(2);
nu_MLE_MM = x;
logL_MM = -fval;
AIC_MM = 2*(k_MM-logL_MM);

results.nu_MM=[m_MLE_MM,K_MLE_MM];
results.logL_MM=logL_MM;
results.AIC_MM=AIC_MM;

fprintf('m_MLE = %4.2f \n',nu_MLE_MM(1));
fprintf('K_MLE = %4.2f \n',nu_MLE_MM(2));
fprintf('Loglikelihood = %4.2f \n',logL_MM);
fprintf('AIC = %4.2f \n',AIC_MM);

% get GE MM model - best fit parameters: m, K
GE_model_MM = Get_GE_models(actual_moi_finer, nu_MLE_MM, Ncells, which_model);
cell_output_MLE_MM = m_MLE_MM*ilist./(K_MLE_MM+ilist); 
% cell_output_MLE_MM(1) = 0;

results.GE_model_MM=GE_model_MM;
results.cell_output_MLE_MM=cell_output_MLE_MM;

%% Hill model:
which_model = 4;

fprintf('-------------------------------------- \n\n');
fprintf('Hill (Saturating) model \n\n');

% m, K, n
x0_Hill = [783.3832, 4.4587, 1]; % Hill function parameters m_init and K_init
k_Hill = length(x0_Hill);
% options = optimset('MaxFunEvals',2e3,'MaxIter',2e3);
% options=[];
[x, fval, exitflag, output] = fminsearch(@(x)ge_inputoutput_models_negloglikelihood(x, data_hpi_18hrs, actual_moi_array, Ncells, GE_noise, which_model), x0_Hill);
m_MLE_Hill = x(1);
K_MLE_Hill = x(2);
n_MLE_Hill = x(3);
nu_MLE_Hill = x;
logL_Hill = -fval;
AIC_Hill = 2*(k_Hill-logL_Hill);

results.nu_Hill=nu_MLE_Hill;
results.logL_Hill=logL_Hill;
results.AIC_Hill=AIC_Hill;

fprintf('m_MLE = %4.2f \n',nu_MLE_Hill(1));
fprintf('K_MLE = %4.2f \n',nu_MLE_Hill(2));
fprintf('n_MLE = %4.2f \n',nu_MLE_Hill(3));
fprintf('Loglikelihood = %4.2f \n',logL_Hill);
fprintf('AIC = %4.2f \n',AIC_Hill);

% get GE MM model - best fit parameters: m, K
GE_model_Hill = Get_GE_models(actual_moi_finer, nu_MLE_Hill, Ncells, which_model);
cell_output_MLE_Hill = m_MLE_MM*ilist.^n_MLE_Hill./(K_MLE_MM+ilist.^n_MLE_Hill); 
% cell_output_MLE_MM(1) = 0;

results.GE_model_Hill=GE_model_Hill;
results.cell_output_MLE_Hill=cell_output_MLE_Hill;

%% now we can get the deltaAIC with respect to MM
deltaAIC_constant = AIC_constant - AIC_MM;
deltaAIC_linear = AIC_linear - AIC_MM;
deltaAIC_MM = AIC_MM - AIC_MM;
deltaAIC_Hill = AIC_Hill - AIC_MM;

deltaAIC = [deltaAIC_constant,deltaAIC_linear,deltaAIC_MM,deltaAIC_Hill];
results.deltaAIC = deltaAIC;

%% print delta-AIC values
fprintf('-------------------------------------- \n\n');

fprintf('Delta AIC_input-independent = %4.2f \n',deltaAIC_constant);
fprintf('Delta AIC_linear = %4.2f \n',deltaAIC_linear);
fprintf('Delta AIC_MM = %4.2f \n\n',deltaAIC_MM);
fprintf('Delta AIC_Hill = %4.2f \n\n',deltaAIC_Hill);


%% plot data and fitting results
figure(1); set(gcf, 'Position',  [100, 800, 1200, 550])

subplot(1,2,1);
plot(actual_moi(:,1), data_hpi_18hrs(:,1), 'ko','MarkerSize',12,'LineWidth',1.5); hold on;
plot(actual_moi(:,2:3), data_hpi_18hrs(:,2:3), 'ko','MarkerSize',12,'LineWidth',1.5); hold on;
h(1) = plot(actual_moi_finer, GE_model_constant, 'Color',colors_rgb_values(1,:),'LineWidth',2); hold on; %plot(actual_moi, GE_model_const, 'm.','MarkerSize',6); 
h(2) = plot(actual_moi_finer, GE_model_linear, 'Color',colors_rgb_values(2,:),'LineWidth',2); hold on; %plot(actual_moi, GE_model_const, 'm.','MarkerSize',6); 
h(3) = plot(actual_moi_finer, GE_model_MM, 'Color',colors_rgb_values(3,:),'LineWidth',2); hold on; %plot(actual_moi, GE_model_const, 'm.','MarkerSize',6); 
h(4) = plot(actual_moi_finer, GE_model_Hill, 'k--','LineWidth',2); hold on; %plot(actual_moi, GE_model_const, 'm.','MarkerSize',6); 
legend(h,{'Input-independent','Linear','Michaelis-Menten','Hill'},'Location','NorthWest');
legend boxoff
xlabel('mean WT MOI'); ylabel('Total Viral Yield, GE/mL')
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 18;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

subplot(1,2,2);
plot(ilist, cell_output_MLE_constant, 'Color',colors_rgb_values(1,:),'LineWidth',2); hold on; 
plot(ilist, cell_output_MLE_constant, '.','MarkerSize',18, 'Color',colors_rgb_values(1,:)); 
plot(ilist, cell_output_MLE_linear, 'Color',colors_rgb_values(2,:),'LineWidth',2); hold on; 
plot(ilist, cell_output_MLE_linear,'.','MarkerSize',18,'Color',colors_rgb_values(2,:)); 
plot(ilist, cell_output_MLE_MM, 'Color',colors_rgb_values(3,:),'LineWidth',2); hold on; 
plot(ilist, cell_output_MLE_MM, '.','MarkerSize',18,'Color',colors_rgb_values(3,:)); 
plot(ilist, cell_output_MLE_Hill, 'k--','LineWidth',2); hold on; 
plot(ilist, cell_output_MLE_Hill, 'k.','MarkerSize',18); 
xlabel('Cellular input (WT Virus)'); ylabel('Cellular output (Total Virus)');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 18;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

%% save file??
if save_file_ans
    
    % save for plotting
    folder_location = '../../Code_plt_ms_figures/results/';
    save(strcat(folder_location,filename),'results','data');
    
    fprintf('File saved:\n'); 
    fprintf(strcat(filename,'\n\n'));
    
    fprintf('Location:\n'); 
    fprintf(strcat(folder_location,'\n\n'));
    
    % save for further analysis
    folder_location = './results/';
    save(strcat(folder_location,filename),'results','data');
    
    fprintf('File saved:\n'); 
    fprintf(strcat(filename,'\n\n'));
    
    fprintf('Location:\n'); 
    fprintf(strcat(folder_location,'\n\n'));
    
    % save for further analysis
    folder_location = '../modelling_TCID50/results/';
    save(strcat(folder_location,filename),'results','data');
    
    fprintf('File saved:\n'); 
    fprintf(strcat(filename,'\n\n'));
    
    fprintf('Location:\n'); 
    fprintf(strcat(folder_location,'\n\n'));
    
    else
    
    fprintf('File not saved.\n');
    
end
