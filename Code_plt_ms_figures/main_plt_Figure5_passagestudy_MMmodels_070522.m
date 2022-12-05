% function void = main_plt_Figure1_passagestudy_2lines_051322(void)

% plot passage study lines 1,2: total particles/TCID50 over 73 passages

%%
clear all; close all; clc;

% save figure?
save_ans_Fig = 0;

figure_name = 'Figure5_passagestudy_2lines_MMmodels_070522';

fprintf('plotting two lines of passage study... \n\n');

color_violet=[169,90,161]/255;
color_lightblue=[133,192,249]/255;
color_orange = [245,121,58]/255;
color_palette = [color_lightblue; color_orange];

frac_spacing = 0.66;
frac_scaling = 0.27;

%%
% load data file (updated 09/10/19)
load('data/Brooke_lines1and2_thru_passage73.mat');
passages = 1:length(params.TCID50); % 73 passages

%%
HAU = params.HAU; % total particles
HAU_LOD = params.HAU_LOD; % 1 - limit of detection (LOD)
HAU_line1 = HAU(1,:);
HAU_line2 = HAU(2,:);

TCID50 = params.TCID50; % 26 - LOD for TCID50
TCID50_LOD= params.TCID50_LOD;
TCID50_line1=TCID50(1,:);
TCID50_line2=TCID50(2,:);

HAU_over_TCID50 = params.HAU./params.TCID50;

ind_HAU_lod_line1 = find(HAU_line1 == HAU_LOD);
ind_HAU_lod_line2 = find(HAU_line2 == HAU_LOD);
HAU_lod_pts_line1 = HAU_line1(ind_HAU_lod_line1);
HAU_lod_pts_line2 = HAU_line2(ind_HAU_lod_line2);

ind_TCID50_lod_line1 = find(TCID50_line1 == TCID50_LOD);
ind_TCID50_lod_line2 = find(TCID50_line2 == TCID50_LOD);
TCID50_lod_pts_line1 = TCID50_line1(ind_TCID50_lod_line1);
TCID50_lod_pts_line2 = TCID50_line2(ind_TCID50_lod_line2);

ind_HAU_TCID50_lod_line1 = [ind_HAU_lod_line1,ind_TCID50_lod_line1];
ind_HAU_TCID50_lod_line2 = [ind_HAU_lod_line2,ind_TCID50_lod_line2];

HAU_over_TCID50_line1 = HAU_line1./TCID50_line1;
HAU_over_TCID50_lod_pts_line1 = HAU_over_TCID50_line1(ind_HAU_TCID50_lod_line1);
HAU_over_TCID50_line2 = HAU_line2./TCID50_line2;
HAU_over_TCID50_lod_pts_line2 = HAU_over_TCID50_line2(ind_HAU_TCID50_lod_line2);

%% now load passage study - MM models
load('results/Results_passagestudy_MMmodels_070522.mat');

HAU_m7 = results.HAU_m7;
HAU_m8 = results.HAU_m8;

TCID50_m7 = results.TCID50_m7;
TCID50_m8 = results.TCID50_m8;



%% now plot things
f1 = figure(1); set(f1, 'Position', [100 500 800 650]);

%% passage study lines 1 and 2
figure(1)
subplot(3,2,1);
semilogy(passages,HAU_LOD*ones(size(passages)),'--k','LineWidth',1); hold on;
p(1)=semilogy(passages, HAU_line1, '.-','Color',[0 0 0],'MarkerSize',14,'LineWidth',2);hold on;
semilogy(ind_HAU_lod_line1, HAU_lod_pts_line2, 'x','Color',[0 0 0],'MarkerSize',8,'LineWidth',2);hold on;
p(2)=semilogy(passages, HAU_line2, '.-','Color',[0.5 0.5 0.5],'MarkerSize',14,'LineWidth',2);hold on;
semilogy(ind_HAU_lod_line2, HAU_lod_pts_line2, 'x','Color',[0.5 0.5 0.5],'MarkerSize',8,'LineWidth',2);hold on;
xlabel('Passage number'); ylabel('HAU');
axis([0 75 5*10^-1 10^3]);
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 16;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

legend(p,{'Line 1','Line 2'},'Location',[0.15, 0.888, .1, .04],'FontSize',10);
legend boxoff;

old_pos = get(f1, 'Position');
set(f1,'Position',[old_pos(1),frac_spacing,old_pos(3),frac_scaling]);
old_pos = get(f1, 'Position');

xticks([0:20:60]);
yticks([10^-1 10^0 10^1 10^2 10^3]);
set(f1,'xticklabel',{[]},'yticklabel',[{' '},{'10^{0}'},{'10^{1}'},{'10^{2}'},{'10^3'}]);

txt = {'A'};
text(0.025,1.05,txt,'Units','normalized','FontSize',14,'FontWeight','bold');

subplot(3,2,3);
semilogy(passages,TCID50_LOD*ones(size(passages)),'--k','MarkerSize',14,'LineWidth',0.5); hold on;
semilogy(passages, TCID50_line1, '.-','Color',[0 0 0],'MarkerSize',14,'LineWidth',2);hold on;
semilogy(ind_TCID50_lod_line1, TCID50_lod_pts_line1, 'x','Color',[0 0 0],'MarkerSize',8,'LineWidth',2);hold on;
semilogy(passages, TCID50_line2, '.-','Color',[0.5 0.5 0.5],'MarkerSize',14,'LineWidth',2);hold on;
semilogy(ind_TCID50_lod_line2, TCID50_lod_pts_line2, 'x','Color',[0.5 0.5 0.5],'MarkerSize',8,'LineWidth',2);hold on;
xlabel('Passage number'); ylabel('TCID_{50}');
axis([0 75 10^0 10^10]);
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 16;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

set(f1,'Position',[old_pos(1), old_pos(2)-frac_scaling, old_pos(3), old_pos(4)])
old_pos = get(f1, 'Position');

xticks([0:20:60]);
yticks([10^0 10^2 10^4 10^6 10^8 10^10]);
set(f1,'xticklabel',{[]},'yticklabel',[{' '},{''},{'10^{4}'},{''},{'10^8'},{' '}]);

txt = {'10^{0}'};
text(-0.1,0.04,txt,'Units','normalized',...
    'FontSize',16,'FontWeight','normal','FontName', 'Times');

subplot(3,2,5);
semilogy(passages, HAU_over_TCID50_line1, '.-','Color',[0 0 0],'MarkerSize',14,'LineWidth',2);hold on;
semilogy(ind_HAU_TCID50_lod_line1, HAU_over_TCID50_lod_pts_line1, 'x','Color',[0 0 0],'MarkerSize',8,'LineWidth',2);hold on;
semilogy(passages, HAU_over_TCID50_line2, '.-','Color',[0.5 0.5 0.5],'MarkerSize',14,'LineWidth',2);hold on;
semilogy(ind_HAU_TCID50_lod_line2, HAU_over_TCID50_lod_pts_line2, 'x','Color',[0.5 0.5 0.5],'MarkerSize',8,'LineWidth',2);hold on;
xlabel('Passage number'); ylabel('HAU/TCID_{50}');
axis([0 75 10^-8 10^1]);
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 16;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

xticks([0:20:60]);
yticks([10^-8 10^-6 10^-4 10^-2 10^0]);
set(f1,'yticklabel',[{'10^{-8}'},{''},{'10^{-4}'},{''},{'10^0'}]);

set(f1,'Position',[old_pos(1), old_pos(2)-frac_scaling, old_pos(3), old_pos(4)])
old_pos = get(f1, 'Position');


%% now plot models
subplot(3,2,2);
q(1) = semilogy(params.passages, HAU_m7*params.HAU_particles_factor, '.-','Color',color_palette(1,:),'MarkerSize',15,'LineWidth',2);hold on;
q(2) = semilogy(params.passages, HAU_m8*params.HAU_particles_factor, '.-','Color',color_palette(2,:),'MarkerSize',15,'LineWidth',2);hold on;
axis([0 75 5*10^5 10^10]);
xlabel('Passage Number'); ylabel('Particles');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 16;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

legend(q,{'Frequency-independent','Frequency-dependent'},'Location',[0.62, 0.888, .1, .04],'FontSize',10);
legend boxoff;

old_pos = get(f1, 'Position');
set(f1,'Position',[old_pos(1),frac_spacing,old_pos(3),frac_scaling]);
old_pos = get(f1, 'Position');

xticks([0:20:60]);
yticks([10^6 10^7 10^8 10^9 10^10]);
set(f1,'xticklabel',{[]},'yticklabel',[{' '},{''},{'10^{8}'},{''},{'10^{10}'}]);

txt = {'10^{6}'};
text(-0.1,0.04,txt,'Units','normalized',...
    'FontSize',16,'FontWeight','normal','FontName', 'Times');

txt = {'B'};
text(0.025,1.05,txt,'Units','normalized','FontSize',14,'FontWeight','bold');


subplot(3,2,4);
semilogy(params.passages, TCID50_m7, '.-','Color',color_palette(1,:),'MarkerSize',15,'LineWidth',2);hold on
semilogy(params.passages, TCID50_m8, '.-','Color',color_palette(2,:),'MarkerSize',15,'LineWidth',2);hold on
axis([0 75 10^0 10^10]);
xlabel('Passage Number'); ylabel('TCID_{50}');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 16;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

set(f1,'Position',[old_pos(1), old_pos(2)-frac_scaling, old_pos(3), old_pos(4)])
old_pos = get(f1, 'Position');

xticks([0:20:60]);
yticks([10^0 10^2 10^4 10^6 10^8 10^10]);
set(f1,'xticklabel',{[]},'yticklabel',[{' '},{''},{'10^{4}'},{''},{'10^{8}'}]);

txt = {'10^{0}'};
text(-0.1,0.04,txt,'Units','normalized',...
    'FontSize',16,'FontWeight','normal','FontName', 'Times');


subplot(3,2,6);
semilogy(params.passages, params.HAU_particles_factor*HAU_m7./TCID50_m7, '.-','Color',color_palette(1,:),'MarkerSize',15,'LineWidth',2);hold on;
semilogy(params.passages, params.HAU_particles_factor*HAU_m8./TCID50_m8, '.-','Color',color_palette(2,:),'MarkerSize',15,'LineWidth',2);hold on;
axis([0 75 10^0 10^9]);
xlabel('Passage Number'); ylabel('Particles/TCID_{50}');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 16;
f1.FontWeight = 'normal';
f1.FontName = 'Times';



set(f1,'Position',[old_pos(1), old_pos(2)-frac_scaling, old_pos(3), old_pos(4)])
old_pos = get(f1, 'Position');

xticks([0:20:60]);
yticks([10^0 10^2 10^4 10^6 10^8]);
set(f1,'yticklabel',[{'10^0'},{''},{'10^{4}'},{''},{'10^{8}'}]);




%
% xlabel('Passage Number'); ylabel('\mbox{HAU}/\mbox{TCID}$_{50}$','Interpreter','Latex');
% axis([0 75 10^-8 10^1]);
% f1=gca;
% f1.LineWidth = 1;
% f1.FontSize = 18;
% f1.FontWeight = 'normal';
% f1.FontName = 'Times';
%
% %%
% yticks([10^-8,10^-6,10^-4,10^-2,10^0]);
% set(f1,'yticklabel',[{'10^{-8}'},{'10^{-6}'},{'10^{-4}'},{'10^{-2}'},{'10^0'}]);
% % legend(g,{'Line 1','Line 2'},'FontSize',18,'Location','Northwest');
% legend(g,{'Line 1','Line 2'},'FontSize',18,'Position', [0.16 0.83 0.12 0.04]);
% % legend(this_p,{legend_char1,legend_char2,legend_char3}, 'Interpreter','Latex','Position', [0.22 0.92 0.098 0.04],'FontSize',9);
% legend boxoff



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