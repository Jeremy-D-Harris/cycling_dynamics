% function void = main_plt_Figure1_passagestudy_2lines_051322(void)

% plot passage study lines 1,2: total particles/TCID50 over 73 passages

%%
clear all; close all; clc;

% save figure?
save_ans_Fig = 0;

figure_name = 'Figure6_passagestudy_plane_2lines_models_080422';

fprintf('plotting two lines of passage study... \n\n');

color_violet=[169,90,161]/255;
color_lightblue=[133,192,249]/255;
color_orange = [245,121,58]/255;
color_palette = [color_lightblue; color_orange];

color_palette_lines12 = [0 0 0; 0.5 0.5 0.5]; % black, gray

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
load('results/Results_passagestudy_allmodels_070722.mat');

HAU_m1 = results.HAU_m1;
HAU_m2 = results.HAU_m2;
HAU_m3 = results.HAU_m3;
HAU_m4 = results.HAU_m4;
HAU_m5 = results.HAU_m5;
HAU_m6 = results.HAU_m6;
HAU_m7 = results.HAU_m7;
HAU_m8 = results.HAU_m8;

TCID50_m1 = results.TCID50_m1;
TCID50_m2 = results.TCID50_m2;
TCID50_m3 = results.TCID50_m3;
TCID50_m4 = results.TCID50_m4;
TCID50_m5 = results.TCID50_m5;
TCID50_m6 = results.TCID50_m6;
TCID50_m7 = results.TCID50_m7;
TCID50_m8 = results.TCID50_m8;


%% limits of detection
ind_HAU_LOD_pts_line2 = find(HAU_line2 == HAU_LOD);
ind_TCID50_LOD_pts_line2 = find(TCID50_line2 == TCID50_LOD);
ind_HAU_TCID50_LOD_pts_line2 = [ind_HAU_LOD_pts_line2,ind_TCID50_LOD_pts_line2];

HAU_LOD_pts_line2 = HAU_line2(ind_HAU_TCID50_LOD_pts_line2);
TCID50_LOD_pts_line2 = TCID50_line2(ind_HAU_TCID50_LOD_pts_line2);

ind_HAU_LOD_pts_line1 = find(HAU_line1 == HAU_LOD);
ind_TCID50_LOD_pts_line1 = find(TCID50_line1 == TCID50_LOD);
ind_HAU_TCID50_LOD_pts_line1 = [ind_HAU_LOD_pts_line1,ind_TCID50_LOD_pts_line1];

HAU_LOD_pts_line1 = HAU_line1(ind_HAU_TCID50_LOD_pts_line1);
TCID50_LOD_pts_line1 = TCID50_line1(ind_HAU_TCID50_LOD_pts_line1);

%% now plot things

f1 = figure(1);
set(gcf, 'Position',  [200, 300, 1550, 400]);

%% panel A: passage study lines 1,2
subplot(1,3,1);

loglog(linspace(1,1e10,length(passages)),HAU_LOD*ones(length(passages)),'--k','MarkerSize',14,'LineWidth',0.5); hold on;
loglog(TCID50_LOD*ones(length(passages)),linspace(1e-1,1e6,length(passages)),'--k','MarkerSize',14,'LineWidth',0.5); hold on;

p(1) = loglog(TCID50_line1, HAU_line1, '-','Color', color_palette_lines12(1,:),'MarkerSize',14,'LineWidth',2);hold on
p(1).Color(4)=0.80;
loglog(TCID50_line1, HAU_line1, 'k.','MarkerSize',14,'LineWidth',2);hold on
loglog(TCID50_LOD_pts_line1, HAU_LOD_pts_line1, 'kx', 'MarkerSize',10,'LineWidth',1.5);hold on

p(2) = loglog(TCID50_line2, HAU_line2, '-','color',color_palette_lines12(2,:),'MarkerSize',14,'LineWidth',2);hold on
p(2).Color(4)=0.80;
loglog(TCID50_line2, HAU_line2, '.','color',color_palette_lines12(2,:),'MarkerSize',14,'LineWidth',2);hold on
loglog(TCID50_LOD_pts_line2, HAU_LOD_pts_line2, 'x', 'MarkerSize',10,'color',color_palette_lines12(2,:),'LineWidth',1.5);hold on


legend(p,{'Line 1', 'Line 2'},'Location',[0.18 0.8 0.05 0.1])
legend boxoff
xlabel('TCID_{50}'); ylabel('HAU');

axis([1 10^10 10^-1 5*10^5]);
% yticks([10^-1 10^0 10^1 10^2 10^3]);
% yticklabels({'10^{-1}', '10^0', '10^1', '10^2', '10^3', '10^4', '10^5'});
ax = gca; % current axes
ax.FontSize = 18;
ax.FontWeight = 'normal';
ax.FontName = 'Times';

txt = {'A'};
text(0.045,1.03,txt,'Units','normalized','FontSize',16,'FontWeight','bold');

%% Panel B: MM models

subplot(1,3,2);
q(1) = loglog(TCID50_m7(100:end), HAU_m7(100:end)*params.HAU_particles_factor, '.-','Color',color_palette(1,:),'MarkerSize',15,'LineWidth',2);hold on
q(2) = loglog(TCID50_m8(100:end), HAU_m8(100:end)*params.HAU_particles_factor, '.-','Color',color_palette(2,:),'MarkerSize',15,'LineWidth',2);hold on
xlabel('TCID_{50}'); ylabel('Particles');
axis([1 10^10 10^6 5e12]);
yticks([10^6 10^7 10^9 10^11]);
yticklabels({ '','10^7', '10^9', '10^{11}'});

% legend(q,{'Line 1', 'MM Frequency-independent', 'MM Frequency-dependent'},'Location',[0.5 0.8 0.1 0.1])
legend(q,{'Frequency-independent', 'Frequency-dependent'},'Location',[0.43 0.8 0.1 0.1]);

% legend(q,{'Saturating frequency-independent', 'Saturating frequency-dependent'},'Location','NorthEast')
% legend(q,{'MM Frequency-independent', 'MM Frequency-dependent'},'Location','NorthEast')
legend boxoff
ax = gca; % current axes
ax.FontSize = 18;
ax.FontWeight = 'normal';
ax.FontName = 'Times';

txt = {'B'};
text(0.045,1.03,txt,'Units','normalized','FontSize',16,'FontWeight','bold');

%% Panel C: models 1-4

subplot(1,3,3);
r(1) = loglog(TCID50_m1(100:end), HAU_m1(100:end)*params.HAU_particles_factor, '.-','Color',color_palette(1,:),'MarkerSize',15,'LineWidth',2);hold on
r(2) = loglog(TCID50_m2(100:end), HAU_m2(100:end)*params.HAU_particles_factor, '.-','Color',color_palette(2,:),'MarkerSize',15,'LineWidth',2);hold on
loglog(TCID50_m3(100:end), HAU_m3(100:end)*params.HAU_particles_factor, '.-','Color',color_palette(1,:),'MarkerSize',15,'LineWidth',2);hold on
loglog(TCID50_m4(100:end), HAU_m4(100:end)*params.HAU_particles_factor, '.-','Color',color_palette(2,:),'MarkerSize',15,'LineWidth',2);hold on
xlabel('TCID_{50}'); ylabel('Particles');
% yticks([10^6 10^7 10^8 10^9 10^10 10^11 10^12]);
% yticklabels({'10^6', '10^7', '10^8', '10^9', '10^{10}', '10^{11}', '10^{12}'});
axis([1 10^10 10^6 5e12]);
% yticks([10^6 10^7 10^9 10^11]);
% yticklabels({ '','10^7', '10^9', '10^{11}'});
% legend(q,{'Line 1', 'MM Frequency-independent', 'MM Frequency-dependent'},'Location',[0.5 0.8 0.1 0.1])

legend(r,{'Frequency-independent', 'Frequency-dependent'},'Location',[0.71 0.8 0.1 0.1]);
legend boxoff
ax = gca; % current axes
ax.FontSize = 18;
ax.FontWeight = 'normal';
ax.FontName = 'Times';

txt = {'C'};
text(0.045,1.03,txt,'Units','normalized','FontSize',16,'FontWeight','bold');


%% place text for labels
text_g1 = {'Input-independent'};
text_g2 = {'Linear'};

text(0.58,0.4,text_g1,'Units','normalized','FontSize',12,'FontWeight','normal','FontName','Times');
text(0.12,0.32,text_g1,'Units','normalized','FontSize',12,'FontWeight','normal','FontName','Times');
text(0.45,0.1,text_g2,'Units','normalized','FontSize',12,'FontWeight','normal','FontName','Times');
text(0.89,0.85,text_g2,'Units','normalized','FontSize',12,'FontWeight','normal','FontName','Times');

yticks([10^6 10^7 10^9 10^11]);
yticklabels({ '','10^7', '10^9', '10^{11}'});

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