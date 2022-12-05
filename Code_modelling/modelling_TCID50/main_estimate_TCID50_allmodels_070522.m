% function void = main_MMmodel_GEpersample_estimatepars_TCID50_newdata_MOI(void)

%% estimate TCID50 vs. DI
% data: using MA primers, assuming intended ratios are met
% To obtain WT & DI MOIs: used DI and WT stock proportion DI for each segment
% see 'Source_Data_3' excel file

% -- updated 07/05/22 --

clear all; close all; clc;

%%
save_file_ans = 0;
% 0: don't save
% 1: save

filename = 'Results_allmodels_TCID50_WTDIratios_070522';

fprintf('Estimating TCID50 vs. WT:DI ratios... \n\n');

% plotting stuff
color_rgb_values = [0 0.4 0.7; 0.9 0.4 0.1];
default_color_rgb = [0    0.4470    0.7410; 0.8500    0.3250    0.0980; 0.9290    0.6940    0.1250; 0.4940    0.1840    0.5560; 0.4660    0.6740    0.1880; 0.3010    0.7450    0.9330; 0.6350    0.0780    0.1840];
Marksize = [35 35 35];
datapoint=['.','.','.'];

%%
% 1-2: Constant models:
% 3-4: Linear models:
% 5-6: Hill models:
% 7-8: Michaelis-Menten (MM) models:

% frequency-independent (odd numbers)
% frequency-dependent (even numbers)

which_model = 1:8;
% which_model = 7;

Ncells = 2e6;
% expected number DI to complement and form WT
eps = 1/24.9589;

% for plotting model fits - GE vs. WT MOI
actual_moi_finer = transpose(linspace(0,10,1000));

%% load fitting results of GE vs. WT MOI
load('./results/FittingResults_GEpersample_070122.mat');
clear data

% 1-2: Constant parameters
constant_pars_fit = results.nu_constant; % nu

% 3-4: Linear parameters
linear_pars_fit = results.nu_linear; % m

% 5-6: Hill parameters
Hill_pars_fit = results.nu_Hill; % m, K, n

% 7-8: Michaelis-Menten parameters
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

for ii=1:length(which_model)
    
    this_model_ind = which_model(ii);
    
    %% 1-2: Constant parameters
    if this_model_ind ==1
        
        %%
        fprintf('-------------------------------------- \n\n');
        fprintf('(1) Input-independent viral yield; \n');
        fprintf('Frequency-independent proportion WT: \n\n');
        
        
        y0_m1 = [0.0201,0.0426]; % phi_init,f_init
        num_pars_m1 = length(y0_m1);
        
        [y, fval_model1, exitflag, output] = fminsearch(@(y)TCID50_inputoutput_allmodels_MOI_negloglikelihood_wDIP(y, newdata_TCID50, TCID50_noise, WT_mois_MA, DI_mois_MA,constant_pars_fit,eps,Ncells,this_model_ind), y0_m1);
        %         phi_fit_m7 = y(1);
        %         f_fit_m7 = y(2);
        pars_fit_m1 = y;
        loglikelihood_model1 = - fval_model1;
        AIC_model1 = 2*(num_pars_m1-loglikelihood_model1);
        
        
        % bulk level: TCID50 vs. DI MOI
        for ii =1:length(DIP_mois_finer)
            
            TCID50_model1(ii) =  Get_TCID50_allmodels_MOI_wDIP(WT_averageMOI_array(ii),DIP_mois_finer(ii),constant_pars_fit,pars_fit_m1,eps,Ncells,this_model_ind); % TCID50 at MOI=1
            
        end
        
        % single-cell level: proportion WT vs. DI MOI
        jlist = 0:100;
        propW_model1 = pars_fit_m1(1)*ones(size(jlist));
        propW_model1(1,1) = 1;
        
        fprintf('phi_MLE = %4.4f \n',pars_fit_m1(1));
        fprintf('f_MLE = %4.4f \n',pars_fit_m1(2));
        fprintf('Loglikelihood = %4.4f \n',loglikelihood_model1);
        %         fprintf('AIC = %4.4f \n',AIC_model1);
        
    elseif this_model_ind ==2
        
        %%
        fprintf('-------------------------------------- \n\n');
        fprintf('(2): Input-independent viral yield; \n');
        fprintf('Frequency-dependent proportion WT: \n\n');
        
        
        y0_m2 = [3.7906,0.0389]; % lambda_init,f_init
        num_pars_m2 = length(y0_m2);
        
        [y, fval_model2, exitflag, output] = fminsearch(@(y)TCID50_inputoutput_allmodels_MOI_negloglikelihood_wDIP(y, newdata_TCID50, TCID50_noise, WT_mois_MA, DI_mois_MA, constant_pars_fit, eps,Ncells,this_model_ind), y0_m2);
        %         lambda_fit = y(1);
        %         f_fit_m8 = y(2);
        pars_fit_m2 = y;
        loglikelihood_model2 = - fval_model2;
        AIC_model2 = 2*(num_pars_m2-loglikelihood_model2);
        %         deltaAIC_models78 = AIC_model7 - AIC_model8
        
        % bulk level: TCID50 vs. DI MOI
        for ii =1:length(DIP_mois_finer)
            
            TCID50_model2(ii) =  Get_TCID50_allmodels_MOI_wDIP(WT_averageMOI_array(ii),DIP_mois_finer(ii),constant_pars_fit,pars_fit_m2,eps,Ncells,this_model_ind); % TCID50 at MOI=1
            
        end
        
        % single-cell level: proportion WT vs. DI MOI
        ilist = [0 1  5 10];
        cnt = 1;
        for ii=ilist
            propW_model2(cnt,:) = (ii*ones(size(jlist))+eps*jlist)./(ii*ones(size(jlist))+eps*jlist+pars_fit_m2(1)*(1-eps)*jlist);
            cnt = cnt+1;
        end
        propW_model2(1,1) = 1;
        
        fprintf('lambda_MLE = %4.4f \n',pars_fit_m2(1));
        fprintf('f_MLE = %4.4f \n',pars_fit_m2(2));
        fprintf('Loglikelihood = %4.4f \n',loglikelihood_model2);
        %         fprintf('AIC = %4.4f \n',AIC_model2);
        
        
        %% 3-4: Linear parameters
    elseif this_model_ind ==3
        
        %%
        fprintf('-------------------------------------- \n\n');
        fprintf('(3) Linear viral yield; \n');
        fprintf('Frequency-independent proportion WT: \n\n');
        
        
        y0_m3 = [0.0101,0.0818]; % phi_init,f_init
        num_pars_m3 = length(y0_m3);
        
        [y, fval_model3, exitflag, output] = fminsearch(@(y)TCID50_inputoutput_allmodels_MOI_negloglikelihood_wDIP(y, newdata_TCID50, TCID50_noise, WT_mois_MA, DI_mois_MA,linear_pars_fit,eps,Ncells,this_model_ind), y0_m3);
        %         phi_fit_m7 = y(1);
        %         f_fit_m7 = y(2);
        pars_fit_m3 = y;
        loglikelihood_model3 = - fval_model3;
        AIC_model3 = 2*(num_pars_m3-loglikelihood_model3);
        
        
        % bulk level: TCID50 vs. DI MOI
        for ii =1:length(DIP_mois_finer)
            
            TCID50_model3(ii) =  Get_TCID50_allmodels_MOI_wDIP(WT_averageMOI_array(ii),DIP_mois_finer(ii),linear_pars_fit,pars_fit_m3,eps,Ncells,this_model_ind); % TCID50 at MOI=1
            
        end
        
        % single-cell level: proportion WT vs. DI MOI
        jlist = 0:100;
        propW_model3 = pars_fit_m3(1)*ones(size(jlist));
        propW_model3(1,1) = 1;
        
        fprintf('phi_MLE = %4.4f \n',pars_fit_m3(1));
        fprintf('f_MLE = %4.4f \n',pars_fit_m3(2));
        fprintf('Loglikelihood = %4.4f \n',loglikelihood_model3);
        %         fprintf('AIC = %4.4f \n',AIC_model3);
        
        
    elseif this_model_ind ==4
        
        %%
        fprintf('-------------------------------------- \n\n');
        fprintf('(4) Linear viral yield; \n');
        fprintf('Frequency-dependent proportion WT: \n\n');
        
        
        y0_m4 = [8.1515 ,0.0769]; % lambda_init,f_init
        num_pars_m4 = length(y0_m4);
        
        [y, fval_model4, exitflag, output] = fminsearch(@(y)TCID50_inputoutput_allmodels_MOI_negloglikelihood_wDIP(y, newdata_TCID50, TCID50_noise, WT_mois_MA, DI_mois_MA, linear_pars_fit, eps,Ncells,this_model_ind), y0_m4);
        %         lambda_fit = y(1);
        %         f_fit_m8 = y(2);
        pars_fit_m4 = y;
        loglikelihood_model4 = - fval_model4;
        AIC_model4 = 2*(num_pars_m4-loglikelihood_model4);
        %         deltaAIC_models78 = AIC_model7 - AIC_model8
        
        % bulk level: TCID50 vs. DI MOI
        for ii =1:length(DIP_mois_finer)
            
            TCID50_model4(ii) =  Get_TCID50_allmodels_MOI_wDIP(WT_averageMOI_array(ii),DIP_mois_finer(ii),linear_pars_fit,pars_fit_m4,eps,Ncells,this_model_ind); % TCID50 at MOI=1
            
        end
        
        % single-cell level: proportion WT vs. DI MOI
        ilist = [0 1  5 10];
        cnt = 1;
        for ii=ilist
            propW_model4(cnt,:) = (ii*ones(size(jlist))+eps*jlist)./(ii*ones(size(jlist))+eps*jlist+pars_fit_m4(1)*(1-eps)*jlist);
            cnt = cnt+1;
        end
        propW_model4(1,1) = 1;
        
        fprintf('lambda_MLE = %4.4f \n',pars_fit_m4(1));
        fprintf('f_MLE = %4.4f \n',pars_fit_m4(2));
        fprintf('Loglikelihood = %4.4f \n',loglikelihood_model4);
        %         fprintf('AIC = %4.4f \n',AIC_model4);
        
        
        %% 5-6: Hill parameters
    elseif this_model_ind ==5
        
        %%
        fprintf('-------------------------------------- \n\n');
        fprintf('(5) Hill viral yield; \n');
        fprintf('Frequency-independent proportion WT: \n\n');
        
        
        y0_m5 = [0.0364,0.0052]; % phi_init,f_init
        num_pars_m5 = length(y0_m5);
        
        [y, fval_model5, exitflag, output] = fminsearch(@(y)TCID50_inputoutput_allmodels_MOI_negloglikelihood_wDIP(y, newdata_TCID50, TCID50_noise, WT_mois_MA, DI_mois_MA,Hill_pars_fit,eps,Ncells,this_model_ind), y0_m5);
        %         phi_fit_m7 = y(1);
        %         f_fit_m7 = y(2);
        pars_fit_m5 = y;
        loglikelihood_model5 = - fval_model5;
        AIC_model5 = 2*(num_pars_m5-loglikelihood_model5);
        
        
        % bulk level: TCID50 vs. DI MOI
        for ii =1:length(DIP_mois_finer)
            
            TCID50_model5(ii) =  Get_TCID50_allmodels_MOI_wDIP(WT_averageMOI_array(ii),DIP_mois_finer(ii),Hill_pars_fit,pars_fit_m5,eps,Ncells,this_model_ind); % TCID50 at MOI=1
            
        end
        
        % single-cell level: proportion WT vs. DI MOI
        jlist = 0:100;
        propW_model5 = pars_fit_m5(1)*ones(size(jlist));
        propW_model5(1,1) = 1;
        
        fprintf('phi_MLE = %4.4f \n',pars_fit_m5(1));
        fprintf('f_MLE = %4.4f \n',pars_fit_m5(2));
        fprintf('Loglikelihood = %4.4f \n',loglikelihood_model5);
        %         fprintf('AIC = %4.2f \n',AIC_model5);
        
    elseif this_model_ind ==6
        
        %%
        fprintf('-------------------------------------- \n\n');
        fprintf('(6) Hill viral yield; \n');
        fprintf('Frequency-dependent proportion WT: \n\n');
        
        
        y0_m6 = [1.801,0.0351]; % lambda_init,f_init
        num_pars_m6 = length(y0_m6);
        
        [y, fval_model6, exitflag, output] = fminsearch(@(y)TCID50_inputoutput_allmodels_MOI_negloglikelihood_wDIP(y, newdata_TCID50, TCID50_noise, WT_mois_MA, DI_mois_MA, Hill_pars_fit, eps,Ncells,this_model_ind), y0_m6);
        %         lambda_fit = y(1);
        %         f_fit_m8 = y(2);
        pars_fit_m6 = y;
        loglikelihood_model6 = - fval_model6;
        AIC_model6 = 2*(num_pars_m6-loglikelihood_model6);
        %         deltaAIC_models78 = AIC_model7 - AIC_model8
        
        % bulk level: TCID50 vs. DI MOI
        for ii =1:length(DIP_mois_finer)
            
            TCID50_model6(ii) =  Get_TCID50_allmodels_MOI_wDIP(WT_averageMOI_array(ii),DIP_mois_finer(ii),Hill_pars_fit,pars_fit_m6,eps,Ncells,this_model_ind); % TCID50 at MOI=1
            
        end
        
        % single-cell level: proportion WT vs. DI MOI
        ilist = [0 1  5 10];
        cnt = 1;
        for ii=ilist
            propW_model6(cnt,:) = (ii*ones(size(jlist))+eps*jlist)./(ii*ones(size(jlist))+eps*jlist+pars_fit_m6(1)*(1-eps)*jlist);
            cnt = cnt+1;
        end
        propW_model6(1,1) = 1;
        
        fprintf('lambda_MLE = %4.4f \n',pars_fit_m6(1));
        fprintf('f_MLE = %4.4f \n',pars_fit_m6(2));
        fprintf('Loglikelihood = %4.4f \n',loglikelihood_model6);
        %         fprintf('AIC = %4.4f \n',AIC_model6);
        
        
        
        %% 7-8: Michaelis-Menten parameters
    elseif this_model_ind ==7
        
        
        %%
        fprintf('-------------------------------------- \n\n');
        fprintf('(7) Michaelis-Menten viral yield; \n');
        fprintf('Frequency-independent proportion WT: \n\n');
        
        
        y0_m7 = [0.011,0.057]; % phi_init,f_init
        num_pars_m7 = length(y0_m7);
        
        [y, fval_model7, exitflag, output] = fminsearch(@(y)TCID50_inputoutput_allmodels_MOI_negloglikelihood_wDIP(y, newdata_TCID50, TCID50_noise, WT_mois_MA, DI_mois_MA,MM_pars_fit,eps,Ncells,this_model_ind), y0_m7);
        %         phi_fit_m7 = y(1);
        %         f_fit_m7 = y(2);
        pars_fit_m7 = y;
        loglikelihood_model7 = - fval_model7;
        AIC_model7 = 2*(num_pars_m7-loglikelihood_model7);
        
        
        % bulk level: TCID50 vs. DI MOI
        for ii =1:length(DIP_mois_finer)
            
            TCID50_model7(ii) =  Get_TCID50_allmodels_MOI_wDIP(WT_averageMOI_array(ii),DIP_mois_finer(ii),MM_pars_fit,pars_fit_m7,eps,Ncells,this_model_ind); % TCID50 at MOI=1
            
        end
        
        % single-cell level: proportion WT vs. DI MOI
        jlist = 0:100;
        propW_model7 = pars_fit_m7(1)*ones(size(jlist));
        propW_model7(1,1) = 1;
        
        fprintf('phi_MLE = %4.4f \n',pars_fit_m7(1));
        fprintf('f_MLE = %4.4f \n',pars_fit_m7(2));
        fprintf('Loglikelihood = %4.4f \n',loglikelihood_model7);
        
    else
        
        %%
        fprintf('-------------------------------------- \n\n');
        fprintf('(8) Michaelis-Menten viral yield; \n');
        fprintf('Frequency-dependent proportion WT: \n\n');
        
        
        y0_m8 = [6.95,0.0536]; % lambda_init,f_init
        num_pars_m8 = length(y0_m8);
        
        [y, fval_model8, exitflag, output] = fminsearch(@(y)TCID50_inputoutput_allmodels_MOI_negloglikelihood_wDIP(y, newdata_TCID50, TCID50_noise, WT_mois_MA, DI_mois_MA, MM_pars_fit, eps,Ncells,this_model_ind), y0_m8);
        %         lambda_fit = y(1);
        %         f_fit_m8 = y(2);
        pars_fit_m8 = y;
        loglikelihood_model8 = - fval_model8;
        AIC_model8 = 2*(num_pars_m8-loglikelihood_model8);
        %         deltaAIC_models78 = AIC_model7 - AIC_model8
        
        % bulk level: TCID50 vs. DI MOI
        for ii =1:length(DIP_mois_finer)
            
            TCID50_model8(ii) =  Get_TCID50_allmodels_MOI_wDIP(WT_averageMOI_array(ii),DIP_mois_finer(ii),MM_pars_fit,pars_fit_m8,eps,Ncells,this_model_ind); % TCID50 at MOI=1
            
        end
        
        % single-cell level: proportion WT vs. DI MOI
        ilist = [0 1  5 10];
        cnt = 1;
        for ii=ilist
            propW_model8(cnt,:) = (ii*ones(size(jlist))+eps*jlist)./(ii*ones(size(jlist))+eps*jlist+pars_fit_m8(1)*(1-eps)*jlist);
            cnt = cnt+1;
        end
        propW_model8(1,1) = 1;
        
        fprintf('lambda_MLE = %4.4f \n',pars_fit_m8(1));
        fprintf('f_MLE = %4.4f \n',pars_fit_m8(2));
        fprintf('Loglikelihood = %4.4f \n',loglikelihood_model8);
        %         fprintf('AIC = %4.4f \n',AIC_model8);
        
        
    end
    
    
end


%% deltaAIC with respect to model 8 (MM+freq-dep)
% deltaAIC_model1 = AIC_model1 - AIC_model8;
% deltaAIC_model2 = AIC_model2 - AIC_model8;
% deltaAIC_model3 = AIC_model3 - AIC_model8;
% deltaAIC_model4 = AIC_model4 - AIC_model8;
% deltaAIC_model5 = AIC_model5 - AIC_model8;
% deltaAIC_model6 = AIC_model6 - AIC_model8;
% deltaAIC_model7 = AIC_model7 - AIC_model8;
% deltaAIC_model8 = AIC_model8 - AIC_model8;

% deltaAIC_allmodels = [deltaAIC_model1,deltaAIC_model2,deltaAIC_model3,deltaAIC_model4,deltaAIC_model5,deltaAIC_model6,deltaAIC_model7,deltaAIC_model8];
% results.deltaAIC = deltaAIC_allmodels;
%
% %% print delta-AIC values
% fprintf('-------------------------------------- \n\n');
%
% fprintf('(1) Delta AIC = %4.4f \n',deltaAIC_allmodels(1));
% fprintf('(2) Delta AIC = %4.4f \n',deltaAIC_allmodels(2));
% fprintf('(3) Delta AIC = %4.4f \n',deltaAIC_allmodels(3));
% fprintf('(4) Delta AIC = %4.4f \n',deltaAIC_allmodels(4));
% fprintf('(5) Delta AIC = %4.4f \n',deltaAIC_allmodels(5));
% fprintf('(6) Delta AIC = %4.4f \n',deltaAIC_allmodels(6));
% fprintf('(7) Delta AIC = %4.4f \n',deltaAIC_allmodels(7));
% fprintf('(8) Delta AIC = %4.4f \n',deltaAIC_allmodels(8));


%% collect data and results before plotting
data_TCID50.WT_mois_MA = WT_mois_MA;
data_TCID50.WT_averageMOI = WT_averageMOI;
data_TCID50.DI_mois_MA = DI_mois_MA;
data_TCID50.DI_mois_MA_shift = DI_mois_MA_shift;
data_TCID50.newdata_TCID50 = newdata_TCID50;

% results of model fitting
results.TCID50_model1=TCID50_model1;
results.TCID50_model2=TCID50_model2;
results.TCID50_model3=TCID50_model3;
results.TCID50_model4=TCID50_model4;
results.TCID50_model5=TCID50_model5;
results.TCID50_model6=TCID50_model6;
results.TCID50_model7=TCID50_model7;
results.TCID50_model8=TCID50_model8;

results.DIP_mois_finer=DIP_mois_finer;

results.propW_model1 = propW_model1;
results.propW_model2 = propW_model2;
results.propW_model3 = propW_model3;
results.propW_model4 = propW_model4;
results.propW_model5 = propW_model5;
results.propW_model6 = propW_model6;
results.propW_model7 = propW_model7;
results.propW_model8 = propW_model8;

results.ilist=ilist;
results.jlist=jlist;

results.pars_fit_m1=pars_fit_m1;
results.pars_fit_m2=pars_fit_m2;
results.pars_fit_m3=pars_fit_m3;
results.pars_fit_m4=pars_fit_m4;
results.pars_fit_m5=pars_fit_m5;
results.pars_fit_m6=pars_fit_m6;
results.pars_fit_m7=pars_fit_m7;
results.pars_fit_m7=pars_fit_m7;
results.pars_fit_m8=pars_fit_m8;

%% plotting...
% WT:DI vs. TCID50
figure(1); set(gcf, 'Position',  [100, 200, 825, 800])



for count=1:4
    
    subplot(4,2,2*count-1);
    
    
    for kk=1:length(DI_mois_MA(:,1))
        for ii=1:3
            if ii == 1
                g(kk) = loglog(max(DI_mois_MA_shift(kk,ii),10^-2), newdata_TCID50(kk,ii),datapoint(ii),'Color', default_color_rgb((kk+1),:),'MarkerSize',Marksize(ii),'LineWidth',1.5); hold on;
            else
                loglog(max(DI_mois_MA_shift(kk,ii),10^-2), newdata_TCID50(kk,ii),datapoint(ii),'Color', default_color_rgb((kk+1),:),'MarkerSize',Marksize(ii),'LineWidth',1.5); hold on;
            end
        end
    end
    
    if count == 1
        g(kk+1) = loglog(DIP_mois_finer, TCID50_model1, 'Color', color_rgb_values(1,:),'MarkerSize',12,'LineWidth',2); hold on;
        g(kk+2) = loglog(DIP_mois_finer, TCID50_model2,'Color',color_rgb_values(2,:),'MarkerSize',12,'LineWidth',2);
        title('Input-independent viral yield');
    elseif count == 2
        g(kk+1) = loglog(DIP_mois_finer, TCID50_model3, 'Color', color_rgb_values(1,:),'MarkerSize',12,'LineWidth',2); hold on;
        g(kk+2) = loglog(DIP_mois_finer, TCID50_model4,'Color',color_rgb_values(2,:),'MarkerSize',12,'LineWidth',2);
        title('Linear-dependent viral yield');
    elseif count == 3
        g(kk+1) = loglog(DIP_mois_finer, TCID50_model5, 'Color', color_rgb_values(1,:),'MarkerSize',12,'LineWidth',2); hold on;
        g(kk+2) = loglog(DIP_mois_finer, TCID50_model6,'Color',color_rgb_values(2,:),'MarkerSize',12,'LineWidth',2);
        title('Hill function viral yield');
    else
        g(kk+1) = loglog(DIP_mois_finer, TCID50_model7, 'Color', color_rgb_values(1,:),'MarkerSize',12,'LineWidth',2); hold on;
        g(kk+2) = loglog(DIP_mois_finer, TCID50_model8,'Color',color_rgb_values(2,:),'MarkerSize',12,'LineWidth',2);
        title('Michaelis-Menten viral yield');
    end
    
    axis([5*10^-3 5*10^2 10^5 10^8]);
    xlabel('DIP MOI'); ylabel('TCID50');
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    %     legend(g,{'WT:DI =  1:0','WT:DI =  1:0.1', 'WT:DI =  1:1', 'WT:DI =  1:10','WT:DI =  1:100','Frequency-independent','Frequency-dependent'},'Location','SouthWest','FontSize',12);
    %     legend boxoff;
    
    %% plot DI MOI (j) vs. proportion WT
    subplot(4,2,2*count);
    
    if count == 1
        
        p(1) = loglog(jlist, propW_model1, '-','Color', color_rgb_values(1,:),'MarkerSize',12,'LineWidth',2); hold on;
        loglog(jlist, propW_model1, '.','Color', color_rgb_values(1,:),'MarkerSize',12,'LineWidth',2); hold on;
        p(2) = loglog(jlist, propW_model2(1,:), 'Color', color_rgb_values(2,:),'MarkerSize',12,'LineWidth',2); hold on;
        loglog(jlist, propW_model2(1,:), '.','Color', color_rgb_values(2,:),'MarkerSize',12,'LineWidth',2); hold on;
        for count=2:length(propW_model2(:,1))
            loglog(jlist, propW_model2(count,:), '-','Color', color_rgb_values(2,:),'MarkerSize',12,'LineWidth',2); hold on;
            loglog(jlist, propW_model2(count,:), '.','Color', color_rgb_values(2,:),'MarkerSize',12,'LineWidth',2); hold on;
        end
        
    elseif count == 2
        
        p(1) = loglog(jlist, propW_model3, '-','Color', color_rgb_values(1,:),'MarkerSize',12,'LineWidth',2); hold on;
        loglog(jlist, propW_model3, '.','Color', color_rgb_values(1,:),'MarkerSize',12,'LineWidth',2); hold on;
        p(2) = loglog(jlist, propW_model4(1,:), 'Color', color_rgb_values(2,:),'MarkerSize',12,'LineWidth',2); hold on;
        loglog(jlist, propW_model4(1,:), '.','Color', color_rgb_values(2,:),'MarkerSize',12,'LineWidth',2); hold on;
        for count=2:length(propW_model4(:,1))
            loglog(jlist, propW_model4(count,:), '-','Color', color_rgb_values(2,:),'MarkerSize',12,'LineWidth',2); hold on;
            loglog(jlist, propW_model4(count,:), '.','Color', color_rgb_values(2,:),'MarkerSize',12,'LineWidth',2); hold on;
        end
        
    elseif count == 3
        
        p(1) = loglog(jlist, propW_model5, '-','Color', color_rgb_values(1,:),'MarkerSize',12,'LineWidth',2); hold on;
        loglog(jlist, propW_model5, '.','Color', color_rgb_values(1,:),'MarkerSize',12,'LineWidth',2); hold on;
        p(2) = loglog(jlist, propW_model6(1,:), 'Color', color_rgb_values(2,:),'MarkerSize',12,'LineWidth',2); hold on;
        loglog(jlist, propW_model6(1,:), '.','Color', color_rgb_values(2,:),'MarkerSize',12,'LineWidth',2); hold on;
        for count=2:length(propW_model6(:,1))
            loglog(jlist, propW_model6(count,:), '-','Color', color_rgb_values(2,:),'MarkerSize',12,'LineWidth',2); hold on;
            loglog(jlist, propW_model6(count,:), '.','Color', color_rgb_values(2,:),'MarkerSize',12,'LineWidth',2); hold on;
        end
        
    else
        p(1) = loglog(jlist, propW_model7, '-','Color', color_rgb_values(1,:),'MarkerSize',12,'LineWidth',2); hold on;
        loglog(jlist, propW_model7, '.','Color', color_rgb_values(1,:),'MarkerSize',12,'LineWidth',2); hold on;
        p(2) = loglog(jlist, propW_model8(1,:), 'Color', color_rgb_values(2,:),'MarkerSize',12,'LineWidth',2); hold on;
        loglog(jlist, propW_model8(1,:), '.','Color', color_rgb_values(2,:),'MarkerSize',12,'LineWidth',2); hold on;
        for count=2:length(propW_model8(:,1))
            loglog(jlist, propW_model8(count,:), '-','Color', color_rgb_values(2,:),'MarkerSize',12,'LineWidth',2); hold on;
            loglog(jlist, propW_model8(count,:), '.','Color', color_rgb_values(2,:),'MarkerSize',12,'LineWidth',2); hold on;
        end
    end
    
    xlabel('DI cellular MOI'); ylabel('Mean proportion WT');
    axis([0 100 10^-3 1]);
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    %     legend(p,{'Frequnecy-independent','Frequency-dependent (i = 0, 1, 5, 10)'},'Location','NorthEast','FontSize',12);
    %     legend boxoff
    
end
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
    folder_location = '../modelling_passagestudy/results/';
    save(strcat(folder_location,filename),'results','data_TCID50');
    
    fprintf('File saved:\n');
    fprintf(strcat(filename,'\n\n'));
    
    fprintf('Location:\n');
    fprintf(strcat(folder_location,'\n\n'));
    
else
    
    fprintf('File not saved.\n');
    
end