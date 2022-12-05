% function void = main_plt_FigureS2_allmodels_passagestudy_plane_070622(void)

% plot passage study lines 1,2: total particles/TCID50 over 73 passages

%%
clear all; close all; clc;

% save figure?
save_ans_Fig = 0;

figure_name = 'FigureS1_allmodelfits_TCID50_080422';

fprintf('plotting fitting results of all model combinations ... \n\n');

% color_violet=[169,90,161]/255;
% color_lightblue=[133,192,249]/255;
% color_orange = [245,121,58]/255;
% color_palette = [color_lightblue; color_orange];
% 
% color_palette_lines12 = [0 0 0; 0.5 0.5 0.5]; % black, gray

% plotting stuff
color_rgb_values = [0 0.4 0.7; 0.9 0.4 0.1];
default_color_rgb = [0    0.4470    0.7410; 0.8500    0.3250    0.0980; 0.9290    0.6940    0.1250; 0.4940    0.1840    0.5560; 0.4660    0.6740    0.1880; 0.3010    0.7450    0.9330; 0.6350    0.0780    0.1840];
Marksize = [35 35 35];
datapoint=['.','.','.'];


frac_spacing = 0.66;
frac_scaling = 0.27;

DIP_mois_finer = 10.^linspace(-3,3,100); %10.^linspace(-2,2,10);

%% now load passage study - MM models
load('results/Results_allmodels_TCID50_WTDIratios_070522.mat');

%% collect data and results before plotting
WT_mois_MA = data_TCID50.WT_mois_MA;
WT_averageMOI = data_TCID50.WT_averageMOI;
DI_mois_MA = data_TCID50.DI_mois_MA;
DI_mois_MA_shift = data_TCID50.DI_mois_MA_shift;
newdata_TCID50 = data_TCID50.newdata_TCID50;

% results of model fitting
TCID50_model1=results.TCID50_model1;
TCID50_model2=results.TCID50_model2;
TCID50_model3=results.TCID50_model3;
TCID50_model4=results.TCID50_model4;
TCID50_model5=results.TCID50_model5;
TCID50_model6=results.TCID50_model6;
TCID50_model7=results.TCID50_model7;
TCID50_model8=results.TCID50_model8;

%% now plot models
% 4 panels: input-independent,linear,Hill,MM

figure(1);
set(gcf, 'Position',  [100, 200, 1250, 300]);


for count=1:4
    
    subplot(1,4,count);
    
    
    for kk=1:length(DI_mois_MA(:,1))
        for ii=1:3
            if ii == 1
                g(kk) = loglog(max(DI_mois_MA_shift(kk,ii),10^-2), newdata_TCID50(kk,ii),datapoint(ii),'Color', default_color_rgb((kk+1),:),'MarkerSize',Marksize(ii),'LineWidth',1.5); hold on;
            else
                loglog(max(DI_mois_MA_shift(kk,ii),10^-2), newdata_TCID50(kk,ii),datapoint(ii),'Color', default_color_rgb((kk+1),:),'MarkerSize',Marksize(ii),'LineWidth',1.5); hold on;
                xticks([10^-2 10^-1 10^0 10^1 10^2]);
                xticklabels({'10^{-2}', '', '10^0', '', '10^2'});
            end
        end
    end
    
    if count == 1
        g(kk+1) = loglog(DIP_mois_finer, TCID50_model1, 'Color', color_rgb_values(1,:),'MarkerSize',12,'LineWidth',2); hold on;
        g(kk+2) = loglog(DIP_mois_finer, TCID50_model2,'Color',color_rgb_values(2,:),'MarkerSize',12,'LineWidth',2);
        title('Input-independent');
        
        legend(g,{'WT:DI =  1:0','WT:DI =  1:0.1', 'WT:DI =  1:1', 'WT:DI =  1:10','WT:DI =  1:100','Frequency-independent','Frequency-dependent'},'Location','SouthWest','FontSize',12);
        legend boxoff;
        
        txt = {'A'};
        text(0.045,1.03,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
        
    elseif count == 2
        loglog(DIP_mois_finer, TCID50_model3, 'Color', color_rgb_values(1,:),'MarkerSize',12,'LineWidth',2); hold on;
        loglog(DIP_mois_finer, TCID50_model4,'Color',color_rgb_values(2,:),'MarkerSize',12,'LineWidth',2);
        title('Linear');
        
        txt = {'B'};
        text(0.045,1.03,txt,'Units','normalized','FontSize',14,'FontWeight','bold');

    elseif count == 3
        loglog(DIP_mois_finer, TCID50_model7, 'Color', color_rgb_values(1,:),'MarkerSize',12,'LineWidth',2); hold on;
        loglog(DIP_mois_finer, TCID50_model8,'Color',color_rgb_values(2,:),'MarkerSize',12,'LineWidth',2);
        title('Michaelis-Menten');
        
        txt = {'C'};
        text(0.045,1.03,txt,'Units','normalized','FontSize',14,'FontWeight','bold');

    else
        loglog(DIP_mois_finer, TCID50_model5, 'Color', color_rgb_values(1,:),'MarkerSize',12,'LineWidth',2); hold on;
        loglog(DIP_mois_finer, TCID50_model6,'Color',color_rgb_values(2,:),'MarkerSize',12,'LineWidth',2);
        title('Hill function');
        
        txt = {'D'};
        text(0.045,1.03,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
        
    end
    
%     axis([5*10^-3 5*10^2 10^5 10^8]);
    axis([10^-3 10^3 10^4 10^8]);
    xlabel('DI mean MOI'); ylabel('TCID_{50}');
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    

    
    
end


% xlabel('DI cellular MOI'); ylabel('Mean proportion WT');
% % axis([0 100 10^-3 1]);
% f1=gca;
% f1.LineWidth = 1;
% f1.FontSize = 14;
% f1.FontWeight = 'normal';
% f1.FontName = 'Times';

%     legend(p,{'Frequnecy-independent','Frequency-dependent (i = 0, 1, 5, 10)'},'Location','NorthEast','FontSize',12);
%     legend boxoff



%% save figure??
if save_ans_Fig
    
    folder_location = './../Figures_ms_all/supp/';
    saveas(f1,strcat(folder_location,figure_name),'epsc');
    
    fprintf('Figure saved:\n');
    fprintf(strcat(figure_name,'\n\n'));
    
    fprintf('Location:\n');
    fprintf(strcat(folder_location,'\n\n'));
    
else
    
    fprintf('Figure not saved.\n');
    
end