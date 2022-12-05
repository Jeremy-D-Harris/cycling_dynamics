% function void = main_estimate_DIPs_GE_050622(void)

%%
clear all; close all; clc;

%%
save_file_ans = 0;
% 0: don't save
% 1: save

filename = 'Results_DIPs_GEpersample_062922';

fprintf('estimating number of DIPs equivalent to WT... \n\n');


% colors for data varying DIP vs. WT ~ 1
setof_rgb_values = [ 0.3867    0.6719    0.7422; 0.0588    0.1255    0.5020;  0.7412    0.7216    0.6784; 0.9333    0.2667    0.1843;   0.6627    0.3529    0.6314];
color_DIonly = [0.3750    0.1016    0.2891];
Marksize = [35 35 35];
datapoint=['.','.','.'];

% which model = Michaelis-Menten
which_model = 3;

% number of cells
Ncells = 2e6;

%% load fitting data
% includes experimenatal data - GE/mL vs. WT MOI - that models were fit to
load('./results/FittingResults_GEpersample_062922.mat');

% GE vs. WT MOI - data that GE models were fit to
data_GE = data.GE_output;
actual_moi = data.actual_moi;
actual_moi_finer = results.actual_moi_finer;

pars_fit = results.nu_MM;
m_fit = pars_fit(1);
K_fit = pars_fit(2);

%% load WT:DIP ratios
% GE & TCID50: three replicates assuming intended WT:DI ratios
% load('./data/PR8_WTvDI_update091219_intendedratios.mat');
load('./data/PR8_WTvDI_update062922.mat');


% TCID50 data
newdata_TCID50 = data_MA.TCID50;

% GE output in the presence of DIPs
newdata_GE_WT_DIextant =  data_MA.GE_output_all;
newdata_GE_WT_DIextant_list = reshape(newdata_GE_WT_DIextant,[1,15]);

% MOIs
WT_mois_MA = data_MA.WT_moi;
DIP_mois_MA = max(data_MA.DI_moi,10^-2);

%% average and range of WT MOI
WT_averageMOI = mean2(WT_mois_MA)
min_WT_mois_MA = min(min(WT_mois_MA))
max_WT_mois_MA = max(max(WT_mois_MA))

% range over DIP:WT ratios
logpart = linspace(-3,3,200);
DIP_WT_ratios_finer = 10.^logpart;
WT_averageMOI_array = WT_averageMOI*ones(size(DIP_WT_ratios_finer));
DIP_mois_finer = WT_averageMOI*DIP_WT_ratios_finer;

%% estimate GE WT or DIPs alone
load('./data/GEpersample_WTorDIP_050622.mat');

WT_1_DI_0 = [1.43E+07 1.18E+07 1.37E+07];
WT_0_DI_10 = [5290928 5508908 4603664]; % this is from the spreadsheet: 'CB_DIPvsWT_data_051318'
viral_yield = [WT_1_DI_0;WT_0_DI_10];


%% estimate number of DIPs that equal WT particle
load('./data/proportion_WTDIsegments.mat');

DI_populations = data_segments.DI_populations; % proportion WT, proportion DIP

% proportion wrt each segment
proportion_WT = transpose(DI_populations(:,1));
proportion_DI = transpose(DI_populations(:,2));

n_virions_list = 1:1000;

% create CMF for probability of at least 1 WT
prev_q_n = prod(1-proportion_DI);
for n_virions = n_virions_list
    
    
    this_q_n = prod(1-proportion_DI.^(n_virions+1));
    prob_at_least_1_WT(n_virions) = this_q_n - prev_q_n;
    prev_q_n = this_q_n;
    
end

% create PMF from CMF
pmf_at_least_1_WT = prob_at_least_1_WT/sum(prob_at_least_1_WT);

DIPs_equivalent_1_WT = sum(n_virions_list.*pmf_at_least_1_WT)

eps_estimate = 1/DIPs_equivalent_1_WT

%% evaluate GE models with Michaelis-Menten formulation
% MLEs
GE_model_MM = results.GE_model_MM;
% increasing WT:DIP ratios
GE_model_MM_wDIP = Get_GE_MM_wDIP(WT_averageMOI,DIP_mois_finer, m_fit, K_fit, eps_estimate, Ncells);

DI_averageMOI_100to1 = mean(data_MA.DI_moi(end,:));
for count=1:length(actual_moi_finer)
    GE_model_MM_wDIP_100to1(1,count) = Get_GE_MM_wDIP(actual_moi_finer(count),DI_averageMOI_100to1, m_fit, K_fit, eps_estimate, Ncells);
end

%% create results structure
results_DIPs.DIPs_equivalent_1_WT = DIPs_equivalent_1_WT;
results_DIPs.eps_estimate = eps_estimate;
results_DIPs.GE_model_MM_wDIP = GE_model_MM_wDIP;
results_DIPs.DI_averageMOI_100to1 = DI_averageMOI_100to1;
results_DIPs.GE_model_MM_wDIP_100to1 = GE_model_MM_wDIP_100to1;

%% now show figures
% probability mass function
figure(1);
bar(n_virions_list, pmf_at_least_1_WT);
axis([0 100 0 0.04]);
xlabel('n DI virions'); ylabel('Probability');
title('Probability of at least 1 WT virus via complementation');

%% panel A: GE vs. WT MOI
figure(2);
set(gcf, 'Position',  [100, 100, 1500, 450]);

subplot(1,3,1)
h(1) = plot(actual_moi(:,1), data_GE(:,1), 'ko','MarkerSize',12,'LineWidth',1.5); hold on;
plot(actual_moi(:,2:3), data_GE(:,2:3), 'ko','MarkerSize',12,'LineWidth',1.5); hold on;

for k=1:length(WT_mois_MA(:,1))
    for i=1:3
        if i == 1
            h(k+1) = plot(WT_mois_MA(k,i), newdata_GE_WT_DIextant(k,i), datapoint(i),'Color',setof_rgb_values(k,:),'MarkerSize',Marksize(i),'LineWidth',1.5); hold on;
        else
            plot(WT_mois_MA(k,i), newdata_GE_WT_DIextant(k,i), datapoint(i),'Color',setof_rgb_values(k,:),'MarkerSize',Marksize(i),'LineWidth',1.5); hold on;
        end
    end
end

h(k+2) = plot(actual_moi_finer, GE_model_MM,'k', 'MarkerSize',12,'LineWidth',2.5); hold on; %,'MarkerSize',12,'LineWidth',1.5); hold on; %plot(effective_moi, GE_model_hill, '.'); hold on;
h(k+3) = plot(actual_moi_finer, GE_model_MM_wDIP_100to1,'k--', 'MarkerSize',12,'LineWidth',2.5); hold on; %,'MarkerSize',12,'LineWidth',1.5); hold on; %plot(effective_moi, GE_model_hill, '.'); hold on;
xlabel('Mean WT MOI'); ylabel('Total viral yield (GE/mL)');
f1 = gca; % current axes
f1.FontSize = 18;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

h_legend1 = 'No DI stock';
h_legend2 = 'WT:DI =  1:0';
h_legend3 = 'WT:DI =  1:0.1';
h_legend4 =  'WT:DI =  1:1';
h_legend5 =  'WT:DI =  1:10';
h_legend6 = 'WT:DI =  1:100';
h_legend7 = 'Saturating model';
h_legend8 = ['Model with DI = ', num2str(DI_averageMOI_100to1,'%2.2f')];
lgnd = legend(h,{h_legend1,h_legend2,h_legend3,h_legend4,h_legend5,h_legend6,h_legend7,h_legend8},'Location',[0.065, 0.66, .25, .25]);
set(lgnd,'FontSize',11);
legend boxoff;

subplot(1,3,1);
text(0.02,0.98,'A','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',22,'FontWeight','bold', 'FontName', 'Arial');

%% panel B: this is to create a bar graph in panel B
% WT:DI = 1:0 and WT:DI = 0:10, y axis is viral yield

% [y_ttest,p_value_ttest] = ttest(WT_1_DI_0,WT_0_DI_10,'Alpha',0.05)

% results_DIPs.y_ttest=y_ttest;
% results_DIPs.p_value_ttest=p_value_ttest;

figure(2);subplot(1,3,2);
for ii = 1:2
    
    h=bar([ii-0.25, ii, ii+0.25],viral_yield(ii,:),0.8); hold on;
    if ii == 1
        set(h,'FaceColor',setof_rgb_values(1,:),'linestyle','none');
    else
        set(h,'FaceColor',color_DIonly(1,:),'linestyle','none');
    end
    %     end
end

xticks([1 2]);
xticklabels({'WT:DI = 1:0','WT:DI = 0:10'});
ylabel('Total viral yield (GE/mL)');
f1 = gca; % current axes
f1.FontSize = 18;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

subplot(1,3,2);
text(0.02,0.98,'B','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',22,'FontWeight','bold', 'FontName', 'Arial')


%% panel C: proportion WT
figure(2);subplot(1,3,3);
h = bar(proportion_WT);
set(h,'FaceColor',color_DIonly(1,:),'linestyle','none');
axis([0 9 0 1]);
xticks([1:8]);
xticklabels({'PB2','PB1','PA','HA','NP','NA','M','NS'});
xlabel('Gene segment'); ylabel('Proportion WT');
f1 = gca; % current axes
f1.FontSize = 18;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

subplot(1,3,3);
text(0.02,0.98,'C','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',22,'FontWeight','bold', 'FontName', 'Arial');



%% save file??
if save_file_ans
    
    % save for plotting
    folder_location = '../../Code_plt_ms_figures/results/';
    save(strcat(folder_location,filename),'results','results_DIPs','data','data_segments','data_MA','data_WTorDIP');
    
    fprintf('File saved:\n');
    fprintf(strcat(filename,'\n\n'));
    
    fprintf('Location:\n');
    fprintf(strcat(folder_location,'\n\n'));
    
else
    
    fprintf('File not saved.\n');
    
end


