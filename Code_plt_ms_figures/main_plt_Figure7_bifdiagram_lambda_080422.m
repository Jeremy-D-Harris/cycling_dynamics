% function void = main_plt_Figure7_bifdiagram_lambda_070622(void)

% bifurcation diagram in lambda

%%
clear all; close all; clc;

% save figure?
save_ans_Fig = 0;

figure_name = 'Figure7_bifdiagram_lambda_080422';

fprintf('plotting bifurcation diagram ... \n\n');

color_violet=[169,90,161]/255;
color_lightblue=[133,192,249]/255;
color_orange = [245,121,58]/255;
color_palette = [color_lightblue; color_orange];

color_palette_lines12 = [0 0 0; 0.65 0.65 0.65]; % black, gray

frac_spacing = 0.66;
frac_scaling = 0.27;


%% now load passage study - MM models
load('results/bifurcationdiagram_lambda_070822.mat');

lambda_vals_collect = [4 6.947 12];


for ii = 1:length(lambda_vals_collect)
    
    ind = find(lambda_vector<lambda_vals_collect(ii));
    ind_lambda_vector(ii) = ind(end);
    lambda_vals_collect_actual(ii) = lambda_vector(ind_lambda_vector(ii));
end

%% now plot models
% 2 panels: bifurcation diagram, trajectories in TCID50-HAU plane

figure(1); set(gcf, 'Position',  [200, 200, 825, 350]);

%% panel A
subplot(1,2,1);
semilogy(lambda_vals_collect_actual(1)*ones(1,10),linspace(10^5,5*10^10,10),'--','Color',color_palette_lines12(2,:),'linewidth',1.5); hold on;
semilogy(lambda_vals_collect_actual(2)*ones(1,10),linspace(10^5,5*10^10,10),'--','Color',color_palette(2,:),'linewidth',1.5); hold on;
semilogy(lambda_vals_collect_actual(3)*ones(1,10),linspace(10^5,5*10^10,10),'--','Color',color_palette_lines12(1,:),'linewidth',1.5); hold on;

semilogy(lambda_vector,max_HAU_varylambda*params.HAU_particles_factor,'k','linewidth',2); hold on;
semilogy(lambda_vector,min_HAU_varylambda*params.HAU_particles_factor,'k','linewidth',2);
xlabel('$\lambda$','interpreter','latex'); 
ylabel('min/max of Particles');
axis([0 15 10^5 5*10^10]);
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 16;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

txt = {'A'};
text(0.045,1.035,txt,'Units','normalized','FontSize',14,'FontWeight','bold');

%% panel B
subplot(1,2,2);
h(1)=loglog(TCID50_varylambda(ind_lambda_vector(1),200:end), HAU_varylambda(ind_lambda_vector(1),200:end)*params.HAU_particles_factor, '.-','Color',color_palette_lines12(2,:),'MarkerSize',15,'LineWidth',2);hold on;
h(2)=loglog(TCID50_varylambda(ind_lambda_vector(2),200:end), HAU_varylambda(ind_lambda_vector(2),200:end)*params.HAU_particles_factor, '.-','Color',color_palette(2,:),'MarkerSize',15,'LineWidth',2); hold on;
h(3)=loglog(TCID50_varylambda(ind_lambda_vector(3),200:end), HAU_varylambda(ind_lambda_vector(3),200:end)*params.HAU_particles_factor, '.-','Color',color_palette_lines12(1,:),'MarkerSize',15,'LineWidth',2); hold on;
axis([10^0 10^10 10^5 5*10^10]);

xlabel('TCID_{50}'); 
ylabel('Particles');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 16;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

txt = {'B'};
text(0.045,1.035,txt,'Units','normalized','FontSize',14,'FontWeight','bold');


h_legend1 = [' $\lambda = $ ', num2str(lambda_vals_collect(1),'%2.2f')];
h_legend2 = [' $\lambda = $ ', num2str(lambda_vals_collect(2),'%2.2f')];
h_legend3 = [' $\lambda = $ ', num2str(lambda_vals_collect(3),'%2.1f')];

lgnd = legend(h,{h_legend1,h_legend2,h_legend3},'interpreter','latex','Location',[0.52, 0.78, .25, .1]);
set(lgnd,'FontSize',12);
legend boxoff;


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