% function void = main_MMmodel_GEpersample_estimatepars_TCID50_newdata_MOI(void)

%% estimate TCID50 vs. DI
% data: using MA primers, assuming intended ratios are met
% To obtain WT & DI MOIs: used DI and WT stock proportion DI for each segment
% see 'Source_Data_3' excel file

% -- updated 07/06/22 --


clear all; close all; clc;

%%
save_file_ans = 0;
% 0: don't save
% 1: save


filename = 'Results_MMmodels_TCID50_WTDIratios_070622';

fprintf('Estimating TCID50 vs. WT:DI ratios... \n\n');

% plotting stuff
color_rgb_values = [0 0.4 0.7; 0.9 0.4 0.1];
default_color_rgb = [0    0.4470    0.7410; 0.8500    0.3250    0.0980; 0.9290    0.6940    0.1250; 0.4940    0.1840    0.5560; 0.4660    0.6740    0.1880; 0.3010    0.7450    0.9330; 0.6350    0.0780    0.1840];
Marksize = [35 35 35];
datapoint=['.','.','.'];

% Michaelis-Menten (MM) models with:
% 7: frequency-independent
% 8: frequency-dependent
MM_models = [7 8];

Ncells = 2e6;
% expected number DI to complement and form WT
eps = 1/24.9589;

% for plotting model fits - GE vs. WT MOI
actual_moi_finer = transpose(linspace(0,10,1000));

%% load fitting results of GE vs. WT MOI
load('./results/FittingResults_GEpersample_070122.mat');
clear data

% Michaelis-Menten parameters
MM_pars_fit = results.nu_MM; % m, K

%% load data - WT 1:DIP X data (MA based)
% get WT and DI MOIs from WT+DI experiments...

% assumes DI and WT stock particles adsorb at their intended ratios
load('./data/PR8_WTvDI_update062922.mat');

% fitting to TCID50 data 
newdata_TCID50 = data_MA.TCID50;

% number of data points
size_data = size(newdata_TCID50);
n_pts = size_data(1)*size_data(2);

% noise is parametrized by a single standard deviation value
% TCID50_noise = mean(std(log10(newdata_TCID50')));

% noise is based on triplicate standard deviation 
TCID50_noise = std(log10(newdata_TCID50'));

newdata_TCID50_vector = reshape(newdata_TCID50,[1, n_pts]);

% WT MOIs (MA)
WT_mois_MA = data_MA.WT_moi;
WT_mois_MA_vector = reshape(WT_mois_MA,[1 n_pts]);

% DIP MOIs (MA)
DI_mois_MA = data_MA.DI_moi;
DI_mois_MA_vector = reshape(DI_mois_MA,[1 n_pts]);

% because the DI MOI and TCID50 values are the same between these
% replicates of WT 1:DI 0
DI_mois_MA_shift = max(DI_mois_MA,10^-2);
DI_mois_MA_shift(1,1) = DI_mois_MA_shift(1,1)-0.001;
DI_mois_MA_shift(1,3) = DI_mois_MA_shift(1,3)+0.001;

% because the DI MOI and TCID50 values are the same between these
% replicates of WT 1:DI 10
DI_mois_MA_shift(4,2) = DI_mois_MA_shift(4,2)-0.3;
DI_mois_MA_shift(4,3) = DI_mois_MA_shift(4,3)+0.3;

WT_averageMOI = mean2(WT_mois_MA); % just to check
WT_maxMOI_MA = max(max(WT_mois_MA));
WT_minMOI_MA = min(min(WT_mois_MA));

DIP_mois_finer = 10.^linspace(-3,3,100); %10.^linspace(-2,2,10);
WT_averageMOI_array = WT_averageMOI*ones(size(DIP_mois_finer));

% legend_WTvDI = 'WT 1:DI X'; %{'WT:DI = 1:0','WT:DI = 1:0.1','WT:DI = 1:1','WT:DI = 1:10'};

%% estimate parameters: 
% (7) frequency-independent model
% (8) frequency-dependent model

for n=1:2
    
    model_ind = MM_models(n);
    
    if model_ind ==7
        
        %%
        fprintf('-------------------------------------- \n\n');
        fprintf('Frequency-independent model: \n\n');
        
        
        y0_m7 = [0.0172679624623060,0.0372244831860066]; % phi_init,f_init
        num_pars_m7 = length(y0_m7);

        [y, fval_model7, exitflag, output] = fminsearch(@(y)TCID50_inputoutput_models78_MOI_negloglikelihood_wDIP(y, newdata_TCID50, TCID50_noise, WT_mois_MA, DI_mois_MA,MM_pars_fit,eps,Ncells,model_ind), y0_m7);
        phi_fit = y(1);
        f_fit_m7 = y(2);
        pars_fit_m7 = [phi_fit,f_fit_m7];
        loglikelihood_model7 = - fval_model7;
        AIC_model7 = 2*(num_pars_m7-loglikelihood_model7);
        
        
        % bulk level: TCID50 vs. DI MOI
        for n =1:length(DIP_mois_finer)
            
            TCID50_model7(n) =  Get_TCID50_models78_MOI_wDIP(WT_averageMOI_array(n),DIP_mois_finer(n),MM_pars_fit,pars_fit_m7,eps,Ncells,model_ind); % TCID50 at MOI=1
            
        end
        
        % single-cell level: proportion WT vs. DI MOI
        jlist = 0:100;
        propW_model7 = phi_fit*ones(size(jlist));
        propW_model7(1,1) = 1;
        
        fprintf('phi_MLE = %4.4f \n',pars_fit_m7(1));
        fprintf('f_MLE = %4.4f \n',pars_fit_m7(2));
        fprintf('Loglikelihood = %4.4f \n',loglikelihood_model7);
%         fprintf('AIC = %4.2f \n',AIC_model7);
        
    else
        
        %%
        fprintf('-------------------------------------- \n\n');
        fprintf('Frequency-dependent model: \n\n');
        
        
        y0_m8 = [5.62365114716726,0.0351211081373544]; % lambda_init,f_init
        num_pars_m8 = length(y0_m8);

        [y, fval_model8, exitflag, output] = fminsearch(@(y)TCID50_inputoutput_models78_MOI_negloglikelihood_wDIP(y, newdata_TCID50, TCID50_noise, WT_mois_MA, DI_mois_MA, MM_pars_fit, eps,Ncells,model_ind), y0_m8);
        lambda_fit = y(1);
        f_fit_m8 = y(2);
        pars_fit_m8 = [lambda_fit f_fit_m8];
        loglikelihood_model8 = - fval_model8;
        AIC_model8 = 2*(num_pars_m8-loglikelihood_model8);
%         deltaAIC_models78 = AIC_model7 - AIC_model8
        
        % bulk level: TCID50 vs. DI MOI
        for n =1:length(DIP_mois_finer)
            
            TCID50_model8(n) =  Get_TCID50_models78_MOI_wDIP(WT_averageMOI_array(n),DIP_mois_finer(n),MM_pars_fit,pars_fit_m8,eps,Ncells,model_ind); % TCID50 at MOI=1
            
        end
        
        % single-cell level: proportion WT vs. DI MOI
        ilist = [0 1  5 10];
        cnt = 1;
        for ii=ilist
            propW_model8(cnt,:) = (ii*ones(size(jlist))+eps*jlist)./(ii*ones(size(jlist))+eps*jlist+lambda_fit*(1-eps)*jlist);
            cnt = cnt+1;
        end
        propW_model8(1,1) = 1;
        
        fprintf('lambda_MLE = %4.4f \n',pars_fit_m8(1));
        fprintf('f_MLE = %4.4f \n',pars_fit_m8(2));
        fprintf('Loglikelihood = %4.4f \n',loglikelihood_model8);
%         fprintf('AIC = %4.4f \n',AIC_model8);
        
        
    end
    
    
end


%% now we can get the deltaAIC with respect to model 8 (MM+freq-dep)
deltaAIC_model7 = AIC_model7 - AIC_model8;
deltaAIC_model8 = AIC_model8 - AIC_model8;

deltaAIC_models78 = [deltaAIC_model7,deltaAIC_model8];
results.deltaAIC = deltaAIC_models78;

%% print delta-AIC values
fprintf('-------------------------------------- \n\n');

fprintf('Delta AIC_frequency-independent = %4.4f \n',deltaAIC_models78(1));
fprintf('Delta AIC_frequency-dependent = %4.4f \n',deltaAIC_models78(2));


%% collect data and results before plotting
data_TCID50.WT_mois_MA = WT_mois_MA;
data_TCID50.WT_averageMOI = WT_averageMOI;
data_TCID50.DI_mois_MA = DI_mois_MA;
data_TCID50.DI_mois_MA_shift = DI_mois_MA_shift;
data_TCID50.newdata_TCID50 = newdata_TCID50;

% results of model fitting
results.TCID50_model7=TCID50_model7;
results.TCID50_model8=TCID50_model8;

results.DIP_mois_finer=DIP_mois_finer;

results.propW_model7 = propW_model7;
results.propW_model8 = propW_model8;

results.ilist=ilist;
results.jlist=jlist;

results.pars_fit_m7=pars_fit_m7;
results.pars_fit_m8=pars_fit_m8;

%% plotting...
% WT:DI vs. TCID50
figure(1); set(gcf, 'Position',  [100, 200, 1150, 400])
subplot(1,2,1);

for kk=1:length(DI_mois_MA(:,1))
    for ii=1:3
        if ii == 1
            g(kk) = loglog(max(DI_mois_MA_shift(kk,ii),10^-2), newdata_TCID50(kk,ii),datapoint(ii),'Color', default_color_rgb((kk+1),:),'MarkerSize',Marksize(ii),'LineWidth',1.5); hold on;
        else
            loglog(max(DI_mois_MA_shift(kk,ii),10^-2), newdata_TCID50(kk,ii),datapoint(ii),'Color', default_color_rgb((kk+1),:),'MarkerSize',Marksize(ii),'LineWidth',1.5); hold on;
        end
    end
end


g(kk+1) = loglog(DIP_mois_finer, TCID50_model7, 'Color', color_rgb_values(1,:),'MarkerSize',12,'LineWidth',2); hold on;
g(kk+2) = loglog(DIP_mois_finer, TCID50_model8,'Color',color_rgb_values(2,:),'MarkerSize',12,'LineWidth',2);
axis([5*10^-3 5*10^2 10^5 10^8]);
xlabel('DIP MOI'); ylabel('TCID50');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 18;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

legend(g,{'WT:DI =  1:0','WT:DI =  1:0.1', 'WT:DI =  1:1', 'WT:DI =  1:10','WT:DI =  1:100','Frequency-independent','Frequency-dependent'},'Location','SouthWest','FontSize',12);
legend boxoff;

%% plot DI MOI (j) vs. proportion WT
figure(1); subplot(1,2,2);
p(1) = loglog(jlist, propW_model7, '-','Color', color_rgb_values(1,:),'MarkerSize',12,'LineWidth',2); hold on;
loglog(jlist, propW_model7, '.','Color', color_rgb_values(1,:),'MarkerSize',12,'LineWidth',2); hold on;
p(2) = loglog(jlist, propW_model8(1,:), 'Color', color_rgb_values(2,:),'MarkerSize',12,'LineWidth',2); hold on;
loglog(jlist, propW_model8(1,:), '.','Color', color_rgb_values(2,:),'MarkerSize',12,'LineWidth',2); hold on;
for count=2:length(propW_model8(:,1))
    loglog(jlist, propW_model8(count,:), '-','Color', color_rgb_values(2,:),'MarkerSize',12,'LineWidth',2); hold on;
    loglog(jlist, propW_model8(count,:), '.','Color', color_rgb_values(2,:),'MarkerSize',12,'LineWidth',2); hold on;
end

xlabel('DI cellular MOI'); ylabel('Mean proportion WT');
axis([0 100 10^-3 1]);
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 18;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

legend(p,{'Frequnecy-independent','Frequency-dependent (i = 0, 1, 5, 10)'},'Location','NorthEast','FontSize',12);
legend boxoff


%% save file??
if save_file_ans
    
    % save for plotting
    folder_location = '../../Code_plt_ms_figures/results/';
    save(strcat(folder_location,filename),'results','data_TCID50');
    
    fprintf('File saved:\n'); 
    fprintf(strcat(filename,'\n\n'));
    
    fprintf('Location:\n'); 
    fprintf(strcat(folder_location,'\n\n'));
    
    % save for further analysis
    folder_location = './results/';
    save(strcat(folder_location,filename),'results','data_TCID50');
    
    fprintf('File saved:\n'); 
    fprintf(strcat(filename,'\n\n'));
    
    fprintf('Location:\n'); 
    fprintf(strcat(folder_location,'\n\n'));
    
    % save for further analysis
    folder_location = '../modelling_passagestudy/results/';
    save(strcat(folder_location,filename),'results','data_TCID50');
    
    fprintf('File saved:\n'); 
    fprintf(strcat(filename,'\n\n'));
    
    fprintf('Location:\n'); 
    fprintf(strcat(folder_location,'\n\n'));
    
    else
    
    fprintf('File not saved.\n');
    
end