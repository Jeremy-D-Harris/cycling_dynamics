% function void = main_plt_Figure2_GEpersample_051322(void)

% plot results from fitting models: WT MOI vs. total viral yield

%%
clear all; close all; clc;

% save figure?
save_ans_Fig = 0;

figure_name = 'Figure2_modelfits_GEpersample_070522';

% 'zesty' color palette: https://venngage.com/blog/color-blind-friendly-palette/
color_violet=[169,90,161]/255;
color_lightblue=[133,192,249]/255;
color_orange = [245,121,58]/255;
color_rgb_values = [color_violet; color_orange; color_lightblue];
% color_rgb_values = [0.6627 0.3529 0.6314; 0.9608  0.4745   0.2275; 0.5216    0.7529    0.9765];

Ncells = 2e6;


%% load fitting data
% includes experimenatal data - GE/mL vs. WT MOI - that models were fit to
load('results/FittingResults_GEpersample_062922.mat');


%%
% GE/mL (MA) after 18 hpi
data_hpi_18hrs = data.GE_output;%data.GE_output;
actual_moi = data.actual_moi;
actual_moi_array = actual_moi;

%% loading results
actual_moi_finer = results.actual_moi_finer;
ilist = 0:15; % this should be saved!!

GE_model_const = results.GE_model_constant;
cell_output_constant = results.cell_output_MLE_constant;

GE_model_linear = results.GE_model_linear; %Get_GE_linear(actual_moi_finer, m_fit_linear, Ncells);
cell_output_linear = results.cell_output_MLE_linear;

GE_model_MM = results.GE_model_MM;
cell_output_MM = results.cell_output_MLE_MM;

%% now plot the results
figure(1); set(gcf, 'Position',  [200, 100, 825, 700]);

%% panel A: MOI vs. viral yield
figure(1); subplot(2,2,1);
h(1) = plot(actual_moi(:,1), data_hpi_18hrs(:,1), 'ko','MarkerSize',10,'LineWidth',1.5); hold on;
plot(actual_moi(:,2:3), data_hpi_18hrs(:,2:3), 'ko','MarkerSize',10,'LineWidth',1.5); hold on;

h(2) = plot(actual_moi_finer, GE_model_const, 'Color',color_rgb_values(1,:),'LineWidth',3); hold on;
h(2).Color(4)=0.85;
h(3) = plot(actual_moi_finer, GE_model_linear, 'Color',color_rgb_values(2,:),'LineWidth',3); hold on;
h(3).Color(4)=0.85;
h(4) = plot(actual_moi_finer, GE_model_MM, 'Color',color_rgb_values(3,:),'LineWidth',3); hold on;
h(4).Color(4)=0.85;
legend boxoff
xlabel('WT Mean MOI'); ylabel('Total Viral Yield (GE/mL)')
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

% set(f1,'yticklabel',[{' '},{'10^{-4}'},{'10^{-2}'},{'10^0'}]);

legend(h,{'Bulk Cell Culture Data','Input-independent','Linear','Saturating'},'FontSize',14,'Location',[0.16 0.82 .15 .075]);

txt = {'A'};
text(0.125,1.04,txt,'Units','normalized','FontSize',14,'FontWeight','bold');


%% panel B: cell MOI vs. viral output
figure(1); subplot(2,2,2);
this_p = plot(ilist, cell_output_constant, 'Color',color_rgb_values(1,:),'LineWidth',3); hold on;
this_p.Color(4)=0.85;
this_p = plot(ilist, cell_output_constant, '.','MarkerSize',18, 'Color',color_rgb_values(1,:)); hold on;
this_p.Color(4)=0.85;
this_p = plot(ilist, cell_output_linear, 'Color',color_rgb_values(2,:),'LineWidth',3); hold on;
this_p.Color(4)=0.85;
this_p = plot(ilist, cell_output_linear,'.','MarkerSize',18,'Color',color_rgb_values(2,:)); hold on;
this_p.Color(4)=0.85;
this_p = plot(ilist, cell_output_MM, 'Color',color_rgb_values(3,:),'LineWidth',3); hold on;
this_p.Color(4)=0.85;
this_p = plot(ilist, cell_output_MM, '.','MarkerSize',18,'Color',color_rgb_values(3,:));
this_p.Color(4)=0.85;
axis([0 15 0 1600]);
xlabel('WT Cellular MOI, $i$','interpreter','latex');  
ylabel('Viral Yield, $\nu(i)$','interpreter','latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

yticks([0:200:1600]);
% set(f1,'yticklabel',[{'0'},{'400'},{'800'},{'1200'},{'1600'}]);
set(f1,'yticklabel',[{'0'},{''},{'400'},{''},{'800'},{''},{'1200'},{''},{'1600'}]);


txt = {'B'};
text(0.025,1.04,txt,'Units','normalized','FontSize',14,'FontWeight','bold');


%% panel C: now load and plot heat map and confidencent intervals
load('results/loglikelihood_M_K_070522.mat');

% bounds for colorbar
lowerbound_logL = -500;
upperbound_logL = max_logL;

% heat map of loglikelihood values
figure(1); subplot(2,2,3);
h = pcolor(K_vector, m_vector,loglikelihood_values); hold on;
% h = pcolor(K_vector, m_vector,max(loglikelihood_values,lowerbound_logL)); hold on;
contour(K_vector, m_vector,loglikelihood_values_binary,'k','LineWidth',1); hold on;
% plot(ci_each_M(:,1),M_values_ci,'k','LineWidth',2); hold on;
% plot(ci_each_M(:,2),M_values_ci,'k','LineWidth',2); hold on;
plot(results.nu_MM(2),results.nu_MM(1),'k.','MarkerSize',20);
set(h, 'EdgeColor', 'none');
xlabel('WT Cellular MOI at Half Maximum Viral Yield, $K$','interpreter','latex'); 
ylabel('Maximum Viral Yield, $m$','interpreter','latex');
cb = colorbar;
set(get(cb,'Title'),'String','Loglikelihood','fontsize',14);
caxis([lowerbound_logL upperbound_logL]);
cb.Ticks = [-500,-400,-300,-200,-100,0];
% cb.Limits = [lowerbound_logL upperbound_logL];
% cb.Ticks = [-500,-400,-300];
cb.TickLabels = {'< -500','','-300','','-100','0'};
% cb.TickLabels = [{'< -500','480','460','440','420','400','380','360','340','320','300'}];
colormap('parula');
axis([2 16 700 2100]);
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

xticks([0:5:15]);
set(f1,'xticklabel',[{'0'},{'5'},{'10'},{'15'}]);
yticks([800:200:2000]);
set(f1,'yticklabel',[{'800'},{' '},{'1200'},{' '},{'1600'},{' '},{'2000'}]);

% set(h, 'DefaultSurfaceEdgeAlpha', 0)

txt = {'C'};
text(0.025,1.04,txt,'Units','normalized','FontSize',14,'FontWeight','bold');

%% now plot MLE input-output relation with 95% fill

% fill in confidence from 95% contour of heat map
for kk = 1:length(locs_positive(:,1))
    [K_matrix,m_matrix] = meshgrid(K_vector,m_vector);
    this_m = m_matrix(locs_positive(kk));
    this_K = K_matrix(locs_positive(kk));
    cell_output_m_K(kk,:) = this_m*ilist./(this_K+ilist);
end

rgb_values_fill = color_rgb_values(3,:)*2/3;

figure(1); subplot(2,2,4)
plot(ilist, cell_output_m_K, 'Color',rgb_values_fill,'LineWidth',3); hold on;
plot(ilist, cell_output_MM, 'Color',color_rgb_values(3,:),'LineWidth',3); hold on;
plot(ilist, cell_output_MM, '.','MarkerSize',18,'Color',color_rgb_values(3,:));
axis([0 15 0 1200])
xlabel('WT Cellular MOI, $i$','interpreter','latex'); 
ylabel('Viral Yield, $\nu(i)$','interpreter','latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

yticks([0:200:1600]);
set(f1,'yticklabel',[{'0'},{' '},{'400'},{' '},{'800'},{' '},{'1200'}]);

txt = {'D'};
text(0.025,1.04,txt,'Units','normalized','FontSize',14,'FontWeight','bold');


%% save figure??
if save_ans_Fig
    
    folder_location = './../Figures_ms_all/main/';
    saveas(f1,strcat(folder_location,figure_name),'epsc');
    
    fprintf('Figure saved:\n');
    fprintf(strcat(figure_name,'\n\n'));
    
    fprintf('Location:\n');
    fprintf(strcat(folder_location,'\n\n'));
    
else
    
    fprintf('Figure not saved.\n');
    
end
