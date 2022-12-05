% function void = main_plt_Figure4_TCID50_WTDIratios_070522(void)

% plot results: DI contribution to WT virus and total viral yield

% -- updated on 07/06/22 --

%%
clear all; close all; clc;


% save figure?
save_ans_Fig = 0;

figure_name = 'Figure4_TCID50_WTDIratios_070622';

fprintf('TCID50 vs. WT:DI ratios... \n\n');

% colors for data varying DIP vs. WT ~ 1
setof_rgb_values = [ 0.3867    0.6719    0.7422; 0.0588    0.1255    0.5020; 0.7412    0.7216    0.6784; 0.9333    0.2667    0.1843; 0.6627    0.3529    0.6314];
color_palette = [0.5216    0.7529    0.9765; 0.9608    0.4745    0.2275]; % deep purple - 0.3750    0.1016    0.289
Marksize = [35 35 35];
datapoint=['.','.','.'];

%% load data & results
load('./results/Results_MMmodels_TCID50_WTDIratios_070522.mat');
newdata_TCID50 = data_TCID50.newdata_TCID50;

% number of data points
size_data = size(newdata_TCID50);
n_pts = size_data(1)*size_data(2);

%% data
WT_mois_MA = data_TCID50.WT_mois_MA;
WT_averageMOI = data_TCID50.WT_averageMOI;
DI_mois_MA = data_TCID50.DI_mois_MA;

DI_mois_MA_shift = data_TCID50.DI_mois_MA_shift;

%% modelling results
TCID50_model7 =  results.TCID50_model7;
TCID50_model8 =  results.TCID50_model8;

DIP_mois_finer = results.DIP_mois_finer;

propW_model7 = results.propW_model7;
propW_model8 = results.propW_model8;

ilist=results.ilist;
jlist=results.jlist;

%% load likelihood heat map and 95% contour
load('results/loglikelihood_lambda_f_070622.mat');

% vector of estimated fraction WT along 95% contour
% f_vector = f_vector_model8;

[f_matrix,lambda_matrix] = meshgrid(f_vector,lambda_vector);
% fill in confidence from 95% contour of heat map
for kk = 1:length(locs_positive(:,1))
    
    this_lambda = lambda_matrix(locs_positive(kk));
    this_f = f_matrix(locs_positive(kk));
    %     these_two_pars_fit_m8 = [this_lambda this_f];
    % fill proportion WT with i = 1, 10
    i=1;
    propW_lambda_f_one(kk,:) = (i*ones(size(jlist))+eps*jlist)./(i*ones(size(jlist))+eps*jlist+this_lambda*(1-eps)*jlist);
    i=10;
    propW_lambda_f_ten(kk,:) = (i*ones(size(jlist))+eps*jlist)./(i*ones(size(jlist))+eps*jlist+this_lambda*(1-eps)*jlist);
    
    
end


% confidence shaded for i = 1
propW_lambda_f_one_minc = min(propW_lambda_f_one);
propW_lambda_f_one_maxc = max(propW_lambda_f_one);
jlist_cat = [jlist, fliplr(jlist)];
inBetweenc_one = [propW_lambda_f_one_minc, fliplr(propW_lambda_f_one_maxc)];

% confidence shaded for i = 10
propW_lambda_f_ten_minc = min(propW_lambda_f_ten);
propW_lambda_f_ten_maxc = max(propW_lambda_f_ten);
inBetweenc_ten = [propW_lambda_f_ten_minc, fliplr(propW_lambda_f_ten_maxc)];

%% redefine variables
TCID50_model7 =  results.TCID50_model7;
TCID50_model8 =  results.TCID50_model8;

DIP_mois_finer = results.DIP_mois_finer;

propW_model7 = results.propW_model7;
propW_model8 = results.propW_model8;

ilist=results.ilist;
jlist=results.jlist;

%% now plot the results
f1 = figure(1); 
% set(gcf, 'Position',  [200, 100, 1000, 800]);
set(gcf, 'Position',  [200, 100, 825, 700]);


subplot(2,2,1);
for kk=1:length(DI_mois_MA(:,1))
    for ii=1:3
        if ii == 1
            %             g(k) = loglog(DIP_WT_ratios(k,i), newdata_TCID50_over_GE(k,i),datapoint(i),'Color', default_color_rgb((k+1),:),'MarkerSize',Marksize(i),'LineWidth',1.5); hold on;
            g(kk) = loglog(DI_mois_MA_shift(kk,ii), newdata_TCID50(kk,ii),datapoint(ii),'Color', setof_rgb_values(kk,:),'MarkerSize',Marksize(ii),'LineWidth',1.5); hold on;
        else
            loglog(DI_mois_MA_shift(kk,ii), newdata_TCID50(kk,ii),datapoint(ii),'Color', setof_rgb_values(kk,:),'MarkerSize',Marksize(ii),'LineWidth',1.5); hold on;
        end
    end
end


this_rgb_value = color_palette(1,:);
figure(1); subplot(2,2,1);
g(kk+1) = loglog(DIP_mois_finer, TCID50_model7, 'Color', this_rgb_value,'MarkerSize',12,'LineWidth',2); hold on;

this_rgb_value = color_palette(2,:);%[0.9 0.4 0.1];
g(kk+2) = loglog(DIP_mois_finer, TCID50_model8,'Color',this_rgb_value,'MarkerSize',12,'LineWidth',2);
xlabel('DI Mean MOI'); 
ylabel('Infectious Virus in Viral Yield (TCID_{50} )');
% xticks([10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3]);
axis([10^-3 10^3 10^4 10^8]);
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';


legend(g,{'WT:DI =  1:0','WT:DI =  1:0.1', 'WT:DI =  1:1', 'WT:DI =  1:10','WT:DI =  1:100','Frequency-independent','Frequency-dependent'},'Location',[0.188 0.63 0.1 0.08]);
legend boxoff;

txt = {'A'};
text(0.025,1.035,txt,'Units','normalized','FontSize',14,'FontWeight','bold');


%% proportion WT vs. Celular DI MOI
figure(1); subplot(2,2,2)
h(1) = loglog(jlist, propW_model7, '.-','Color', color_palette(1,:),'MarkerSize',16,'LineWidth',2); hold on;
h(2) = loglog(jlist, propW_model8(1,:),'.-', 'Color', color_palette(2,:),'MarkerSize',16,'LineWidth',2); hold on;

% plot for i=0,1,5,10
for kk=2:length(propW_model8(:,1))
    
    loglog(jlist, propW_model8(kk,:), '.-','Color', color_palette(2,:),'MarkerSize',16,'LineWidth',2); hold on;
    
end
legend(h,{'Frequency-independent','Frequency-dependent'},'Location',[0.62 0.59 0.13 0.05]);
legend boxoff
xlabel('DI Cellular MOI, $j$','interpreter','latex'); 
ylabel('Proportion WT in Viral Yield, $\phi_W(i,j)$','interpreter','latex');
axis([1 100 10^-3 1])
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

txt = {'B'};
text(0.025,1.035,txt,'Units','normalized','FontSize',14,'FontWeight','bold');

%%
subplot(2,2,2)
txt = {'$i = 0$'};
text(0.05,0.22,txt,'interpreter','latex','Units','normalized','FontSize',12,'FontWeight','normal','FontName','Times');
txt = {'$i = 1$'};
text(0.05,0.58,txt,'interpreter','latex','Units','normalized','FontSize',12,'FontWeight','normal','FontName','Times');
txt = {'$i = 5$'};
text(0.05,0.78,txt,'interpreter','latex','Units','normalized','FontSize',12,'FontWeight','normal','FontName','Times');
txt = {'$i = 10$'};
text(0.05,0.94,txt,'interpreter','latex','Units','normalized','FontSize',12,'FontWeight','normal','FontName','Times');


%% panel C: loglikelihood heat map with contour plot
lowerbound_logL = -100;
upperbound_logL = max_logL;

figure(1); subplot(2,2,3);
h = pcolor(f_vector,lambda_vector,loglikelihood_values_model8); hold on;
M=contour(f_vector, lambda_vector,loglikelihood_values_binary,1,'k','LineWidth',1.5); hold on;
% plot(f_max_m8,lambda_max_m8,'k.','MarkerSize',20);
plot(results.pars_fit_m8(2),results.pars_fit_m8(1),'k.','MarkerSize',20);
axis([f_vector(1) f_vector(end) lambda_vector(1) lambda_vector(end)])
set(h, 'EdgeColor', 'none');
xlabel('Fraction of Infectious WT Virus, $f$','interpreter','latex'); 
ylabel({'Fitness of DI relative to WT Virus, $\lambda$'},'interpreter','latex');
cb = colorbar;
set(get(cb,'Title'),'String','Loglikelihood','fontsize',12);
caxis([lowerbound_logL upperbound_logL]);
cb.Ticks = [-100,-80,-60,-40,-20,0];
cb.TickLabels = {'< -100','','-60','','-20','0'};
colormap('parula')
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

%% confidence intervals for f, lambda
f_ub_95percent = max(M(1,2:end))
f_lb_95percent = min(M(1,2:end))

lambda_ub_95percent = max(M(2,2:end))
lambda_lb_95percent = min(M(2,2:end))

%%
txt = {'C'};
text(0.025,1.035,txt,'Units','normalized','FontSize',14,'FontWeight','bold');

%% panel D: frequency-dependent: i=1,10 with 95% intervals
this_rgb_value = color_palette(2,:);%[0.9 0.4 0.1];
rgb_values_fill = this_rgb_value*2/3;%[1 0.7 0.5];

figure(1); subplot(2,2,4);
loglog(jlist, propW_lambda_f_one, 'Color',rgb_values_fill,'LineWidth',2); hold on;
loglog(jlist, propW_lambda_f_ten, 'Color',rgb_values_fill,'LineWidth',2); hold on;
p = loglog(jlist, propW_model8(2,:), '.-', 'Color',this_rgb_value,'MarkerSize',16,'LineWidth',2); hold on;
loglog(jlist, propW_model8(4,:), '.-', 'Color',this_rgb_value,'MarkerSize',16,'LineWidth',2); hold on;
% loglog(jlist,propW_lambda_f_ten_maxc,jlist,propW_lambda_f_ten_minc, 'Color',this_rgb_value,'MarkerSize',16,'LineWidth',1); hold on;
% loglog(jlist,propW_lambda_f_one_maxc,jlist,propW_lambda_f_one_minc, 'Color',this_rgb_value,'MarkerSize',16,'LineWidth',1); hold on;
axis([1 100 10^-3 1]);
xlabel('DI Cellular MOI, $j$','interpreter','latex'); 
ylabel('Proportion WT in Viral Yield, $\phi_W(i,j)$','interpreter','latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';


% legend(p,{'Frequency-dependent, $i = 1, 10$'},'interpreter','latex','Location','Southeast');
legend(p,{'Frequency-dependent ($i = 1,\, 10$) '},'interpreter','latex','Location',[0.675 0.11 0.1 0.05]);
legend boxoff;

txt = {'D'};
text(0.025,1.035,txt,'Units','normalized','FontSize',14,'FontWeight','bold');


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




