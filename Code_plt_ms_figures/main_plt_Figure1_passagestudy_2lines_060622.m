% function void = main_plt_Figure1_passagestudy_2lines_051322(void)

% plot passage study lines 1,2: total particles/TCID50 over 73 passages

%%
clear all; close all; clc;

% save figure?
save_ans_Fig = 0;

figure_name = 'Figure1_passagestudy_lines1and2_060622_gray';

fprintf('plotting two lines of passage study... \n\n');


color_palette = [0 0 0; 0.5 0.5 0.5]; % black, gray
% or use
% same as in figure 3
% violet_rgb = [0.3750    0.1016    0.2891]; 
% teal_rgb = [0.3867    0.6719    0.7422];
% color_palette = [violet_rgb; teal_rgb];

%%
% load data file (updated 09/10/19)
load('data/Brooke_lines1and2_thru_passage73.mat');
passages = 1:length(params.TCID50); % 73 passages

%%
HAU = params.HAU; % total particles
HAU_LOD = params.HAU_LOD; % 1 - limit of detection (LOD)
TCID50 = params.TCID50; % 26 - LOD for TCID50
TCID50_LOD= params.TCID50_LOD;
HAU_over_TCID50 = params.HAU./params.TCID50;

ind_HAU_LOD_line1 = find(HAU(1,:) == HAU_LOD);
ind_HAU_LOD_line2 = find(HAU(2,:) == HAU_LOD);

ind_TCID50_LOD_line1 = find(TCID50(1,:) == TCID50_LOD);
ind_TCID50_LOD_line2 = find(TCID50(2,:) == TCID50_LOD);

ind_HAU_TCID50_LOD_line1 = [ind_HAU_LOD_line1,ind_TCID50_LOD_line1];
ind_HAU_TCID50_LOD_line2 = [ind_HAU_LOD_line2,ind_TCID50_LOD_line2];

%% now plot the lines
figure(1); set(gcf, 'Position',  [500, 800, 875, 600]);
% subplot(1,2,1);

% loglog(ind_HAU_TCID50_LOD_line1,HAU_over_TCID50(1,ind_HAU_TCID50_LOD_line1),'r.','MarkerSize',15,'LineWidth',2); hold on;

g(1) = semilogy(passages,HAU_over_TCID50(1,:),'Color',color_palette(1,:),'LineWidth',2.5); hold on;
g(1).Color(4) = 0.85;
semilogy(passages,HAU_over_TCID50(1,:),'.','Color',color_palette(1,:),'MarkerSize',16); hold on;

g(2) = semilogy(passages,HAU_over_TCID50(2,:),'Color',color_palette(2,:),'LineWidth',2.5); hold on;
g(2).Color(4) = 0.85;
semilogy(passages,HAU_over_TCID50(2,:),'.','Color',color_palette(2,:),'MarkerSize',16); hold on;

% loglog(ind_HAU_TCID50_LOD_line2,HAU_over_TCID50(2,ind_HAU_TCID50_LOD_line2),'r.','MarkerSize',15,'LineWidth',2); hold on;

%%
xlabel('Passage Number'); ylabel('\mbox{HAU}/\mbox{TCID}$_{50}$','Interpreter','Latex');
axis([0 75 10^-8 10^1]);
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 18;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

%%
yticks([10^-8,10^-6,10^-4,10^-2,10^0]);
set(f1,'yticklabel',[{'10^{-8}'},{'10^{-6}'},{'10^{-4}'},{'10^{-2}'},{'10^0'}]);
% legend(g,{'Line 1','Line 2'},'FontSize',18,'Location','Northwest');
legend(g,{'Line 1','Line 2'},'FontSize',18,'Position', [0.16 0.83 0.12 0.04]);
% legend(this_p,{legend_char1,legend_char2,legend_char3}, 'Interpreter','Latex','Position', [0.22 0.92 0.098 0.04],'FontSize',9);
legend boxoff



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