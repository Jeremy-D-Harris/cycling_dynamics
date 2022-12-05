% function void = main_plt_FigureS2_allmodels_passagestudy_plane_070622(void)

% plot passage study lines 1,2: total particles/TCID50 over 73 passages

%%
clear all; close all; clc;

% save figure?
save_ans_Fig = 0;

figure_name = 'FigureS2_allmodels_passagestudy_plane_080422';

fprintf('plotting parametrized models of passage study... \n\n');

color_violet=[169,90,161]/255;
color_lightblue=[133,192,249]/255;
color_orange = [245,121,58]/255;
color_palette = [color_lightblue; color_orange];

color_palette_lines12 = [0 0 0; 0.5 0.5 0.5]; % black, gray

frac_spacing = 0.66;
frac_scaling = 0.27;


%% now load passage study - MM models
load('results/Results_passagestudy_allmodels_070422.mat');

HAU_m1 = results.HAU_m1(1,100:end);
HAU_m2 = results.HAU_m2(1,100:end);
HAU_m3 = results.HAU_m3(1,100:end);
HAU_m4 = results.HAU_m4(1,100:end);
HAU_m5 = results.HAU_m5(1,100:end);
HAU_m6 = results.HAU_m6(1,100:end);
HAU_m7 = results.HAU_m7(1,100:end);
HAU_m8 = results.HAU_m8(1,100:end);

TCID50_m1 = results.TCID50_m1(1,100:end);
TCID50_m2 = results.TCID50_m2(1,100:end);
TCID50_m3 = results.TCID50_m3(1,100:end);
TCID50_m4 = results.TCID50_m4(1,100:end);
TCID50_m5 = results.TCID50_m5(1,100:end);
TCID50_m6 = results.TCID50_m6(1,100:end);
TCID50_m7 = results.TCID50_m7(1,100:end);
TCID50_m8 = results.TCID50_m8(1,100:end);

%% now plot models
% 4 panels: input-independent,linear,Hill,MM

figure(1);
set(gcf, 'Position',  [100, 200, 1250, 300]);

%% input-independent
figure(1); subplot(1,4,1);
p(1)=loglog(TCID50_m1, HAU_m1, '.-','Color',color_palette(1,:),'MarkerSize',15,'LineWidth',2); hold on;
p(2)=loglog(TCID50_m2, HAU_m2, '.-','Color',color_palette(2,:),'MarkerSize',15,'LineWidth',2); hold on;

xlabel('TCID_{50}'); ylabel('HAU');
title('Input-independent');
axis([10^0 10^10 10^-2 10^5]);
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

legend(p,{'Frequency-independent','Frequency-dependent'},'Location',[0.15 0.79 0.13 0.05]);
legend boxoff

txt = {'A'};
text(0.045,1.03,txt,'Units','normalized','FontSize',14,'FontWeight','bold');

%% linear
figure(1); subplot(1,4,2);
loglog(TCID50_m3, HAU_m3, '.-','Color',color_palette(1,:),'MarkerSize',15,'LineWidth',2); hold on;
loglog(TCID50_m4, HAU_m4, '.-','Color',color_palette(2,:),'MarkerSize',15,'LineWidth',2); hold on;

xlabel('TCID_{50}'); ylabel('HAU');
title('Linear');
axis([10^0 10^10 10^-2 10^5]);
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

txt = {'B'};
text(0.045,1.03,txt,'Units','normalized','FontSize',14,'FontWeight','bold');

%% MM
figure(1); subplot(1,4,3);
loglog(TCID50_m7, HAU_m7, '.-','Color',color_palette(1,:),'MarkerSize',15,'LineWidth',2); hold on;
loglog(TCID50_m8, HAU_m8, '.-','Color',color_palette(2,:),'MarkerSize',15,'LineWidth',2); hold on;

xlabel('TCID_{50}'); ylabel('HAU');
title('Michaelis-Menten');
axis([10^0 10^10 10^-2 10^5]);
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

txt = {'C'};
text(0.045,1.03,txt,'Units','normalized','FontSize',14,'FontWeight','bold');

%% Hill
figure(1); subplot(1,4,4);
loglog(TCID50_m5, HAU_m5, '.-','Color',color_palette(1,:),'MarkerSize',15,'LineWidth',2); hold on;
loglog(TCID50_m6, HAU_m6, '.-','Color',color_palette(2,:),'MarkerSize',15,'LineWidth',2); hold on;

xlabel('TCID_{50}'); ylabel('HAU');
title('Hill function');
axis([10^0 10^10 10^-2 10^5]);
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

txt = {'D'};
text(0.045,1.03,txt,'Units','normalized','FontSize',14,'FontWeight','bold');

% 
% %% now plot things
% 
% f1 = figure(1);
% set(gcf, 'Position',  [200, 300, 1550, 400]);
% 
% %% panel A: passage study lines 1,2
% subplot(1,3,1);
% 
% loglog(linspace(1,1e10,length(passages)),HAU_LOD*ones(length(passages)),'--k','MarkerSize',14,'LineWidth',0.5); hold on;
% loglog(TCID50_LOD*ones(length(passages)),linspace(1e-1,1e6,length(passages)),'--k','MarkerSize',14,'LineWidth',0.5); hold on;
% 
% p(1) = loglog(TCID50_line1, HAU_line1, '-','Color', color_palette_lines12(1,:),'MarkerSize',14,'LineWidth',2);hold on
% p(1).Color(4)=0.80;
% loglog(TCID50_line1, HAU_line1, 'k.','MarkerSize',14,'LineWidth',2);hold on
% loglog(TCID50_LOD_pts_line1, HAU_LOD_pts_line1, 'kx', 'MarkerSize',10,'LineWidth',1.5);hold on
% 
% p(2) = loglog(TCID50_line2, HAU_line2, '-','color',color_palette_lines12(2,:),'MarkerSize',14,'LineWidth',2);hold on
% p(2).Color(4)=0.80;
% loglog(TCID50_line2, HAU_line2, '.','color',color_palette_lines12(2,:),'MarkerSize',14,'LineWidth',2);hold on
% loglog(TCID50_LOD_pts_line2, HAU_LOD_pts_line2, 'x', 'MarkerSize',10,'color',color_palette_lines12(2,:),'LineWidth',1.5);hold on
% 
% 
% legend(p,{'Line 1', 'Line 2'},'Location',[0.18 0.8 0.05 0.1])
% legend boxoff
% xlabel('TCID_{50}'); ylabel('HAU');
% 
% axis([1 10^10 10^-1 5*10^5]);
% % yticks([10^-1 10^0 10^1 10^2 10^3]);
% % yticklabels({'10^{-1}', '10^0', '10^1', '10^2', '10^3', '10^4', '10^5'});
% ax = gca; % current axes
% ax.FontSize = 16;
% ax.FontWeight = 'normal';
% ax.FontName = 'Times';
% 
% txt = {'A'};
% text(0.045,1.03,txt,'Units','normalized','FontSize',16,'FontWeight','bold');
% 
% %% Panel B: MM models
% 
% subplot(1,3,2);
% q(1) = loglog(TCID50_m7(100:end), HAU_m7(100:end)*params.HAU_particles_factor, '.-','Color',color_palette(1,:),'MarkerSize',15,'LineWidth',2);hold on
% q(2) = loglog(TCID50_m8(100:end), HAU_m8(100:end)*params.HAU_particles_factor, '.-','Color',color_palette(2,:),'MarkerSize',15,'LineWidth',2);hold on
% xlabel('TCID_{50}'); ylabel('HAU-particles');
% axis([1 10^10 10^6 5e12]);
% yticks([10^6 10^7 10^9 10^11]);
% yticklabels({ '','10^7', '10^9', '10^{11}'});
% 
% % legend(q,{'Line 1', 'MM Frequency-independent', 'MM Frequency-dependent'},'Location',[0.5 0.8 0.1 0.1])
% legend(q,{'Frequency-independent', 'Frequency-dependent'},'Location',[0.43 0.8 0.1 0.1]);
% 
% % legend(q,{'Saturating frequency-independent', 'Saturating frequency-dependent'},'Location','NorthEast')
% % legend(q,{'MM Frequency-independent', 'MM Frequency-dependent'},'Location','NorthEast')
% legend boxoff
% ax = gca; % current axes
% ax.FontSize = 16;
% ax.FontWeight = 'normal';
% ax.FontName = 'Times';
% 
% txt = {'B'};
% text(0.045,1.03,txt,'Units','normalized','FontSize',16,'FontWeight','bold');
% 
% %% Panel C: models 1-4
% 
% subplot(1,3,3);
% r(1) = loglog(TCID50_m1(100:end), HAU_m1(100:end)*params.HAU_particles_factor, '.-','Color',color_palette(1,:),'MarkerSize',15,'LineWidth',2);hold on
% r(2) = loglog(TCID50_m2(100:end), HAU_m2(100:end)*params.HAU_particles_factor, '.-','Color',color_palette(2,:),'MarkerSize',15,'LineWidth',2);hold on
% loglog(TCID50_m3(100:end), HAU_m3(100:end)*params.HAU_particles_factor, '.-','Color',color_palette(1,:),'MarkerSize',15,'LineWidth',2);hold on
% loglog(TCID50_m4(100:end), HAU_m4(100:end)*params.HAU_particles_factor, '.-','Color',color_palette(2,:),'MarkerSize',15,'LineWidth',2);hold on
% xlabel('TCID_{50}'); ylabel('HAU-particles');
% % yticks([10^6 10^7 10^8 10^9 10^10 10^11 10^12]);
% % yticklabels({'10^6', '10^7', '10^8', '10^9', '10^{10}', '10^{11}', '10^{12}'});
% axis([1 10^10 10^6 5e12]);
% % yticks([10^6 10^7 10^9 10^11]);
% % yticklabels({ '','10^7', '10^9', '10^{11}'});
% % legend(q,{'Line 1', 'MM Frequency-independent', 'MM Frequency-dependent'},'Location',[0.5 0.8 0.1 0.1])
% 
% legend(r,{'Frequency-independent', 'Frequency-dependent'},'Location',[0.71 0.8 0.1 0.1]);
% legend boxoff
% ax = gca; % current axes
% ax.FontSize = 16;
% ax.FontWeight = 'normal';
% ax.FontName = 'Times';
% 
% txt = {'C'};
% text(0.045,1.03,txt,'Units','normalized','FontSize',16,'FontWeight','bold');
% 
% 
% %% place text for labels
% text_g1 = {'Input-independent'};
% text_g2 = {'Linear'};
% 
% text(0.58,0.4,text_g1,'Units','normalized','FontSize',12,'FontWeight','normal','FontName','Times');
% text(0.12,0.32,text_g1,'Units','normalized','FontSize',12,'FontWeight','normal','FontName','Times');
% text(0.45,0.1,text_g2,'Units','normalized','FontSize',12,'FontWeight','normal','FontName','Times');
% text(0.89,0.85,text_g2,'Units','normalized','FontSize',12,'FontWeight','normal','FontName','Times');
% 
% yticks([10^6 10^7 10^9 10^11]);
% yticklabels({ '','10^7', '10^9', '10^{11}'});

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