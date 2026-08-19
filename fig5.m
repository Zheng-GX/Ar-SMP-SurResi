close all;clear;clc
% 
set(0,'defaultAxesFontName', 'Times new roman','defaultAxesFontsize',17);
set(0,'defaultTextInterpreter','latex','defaultLegendInterpreter','latex'); 
set(0,'defaultfigurecolor','w','defaultAxeslinewidth',1)
set(0,'defaultfigureposition',[500,300,450,320]);
c1=[91,181,172]/255;
c2=[216,179,101]/255;
c3=[155,155,155]/255;
%%
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% xd=50
load("a6_FPT_20260609_50.mat")
%%%%
sur(1)=1-sum(isnan(FPT))/n;
sur1(1)=1-sum(isnan(FPT1))/n;
lins=':';
%%%%
figure('Position',[732,312,489,395])

ax1=subplot(2,3,1);

% set(gca,'Position',[0.151,0.59,0.262,0.341])
hold on;
%%%%
FPT_clean = FPT(~isnan(FPT));
FPT_clean1 = FPT1(~isnan(FPT1));
%%%%%%%%%% SMP
[f_log, xi_log] = ksdensity(log(FPT_clean), 'NumPoints', 1000, 'Support', 'positive','Bandwidth','plug-in');
% 
xi_orig = exp(xi_log);                      
f_orig = f_log ./ xi_orig;                   
%%%%%%%%%% MP
[f_log1, xi_log1] = ksdensity(log(FPT_clean1), 'NumPoints', 1000, 'Support', 'positive','Bandwidth','plug-in');
% 
xi_orig1 = exp(xi_log1);                       
f_orig1 = f_log1 ./ xi_orig1;                   
%
plot(xi_orig, f_orig, '-', 'LineWidth',2, 'color',c1,'DisplayName','SMP','MarkerSize',3);
plot(xi_orig1, f_orig1, lins, 'LineWidth',2, 'color',c2,'DisplayName','MP');
% fill([xi_orig, fliplr(xi_orig)], [f_orig, zeros(size(f_orig))],color,...
%         'FaceAlpha', 0.9,'EdgeColor', 'none','HandleVisibility', 'off'); 
legend('Box','off','Position',[0.22063773674997,0.658574001021107,0.178255129679168,0.113468357858778])
% text(20,0.02,['$\left \langle \rm FPT \right \rangle =$' num2str(round(mean_FPT)),' s'],'Interpreter','latex','FontSize',15)
% xlabel('FPT (s)')
ylabel('Probability density (1/s)','Position',[-60.041818354676025,-0.004537020273231,-1])
xlim([5 200])
ylim([0 0.035])
% yticks([0 0.005 0.01])
% yticks([0 0.01 0.02 0.03 ])
xticks([50 100 150])
yticks(0:0.01:0.03)
title('$x_d=50\ \mu \rm m$','FontSize',14)
box on
%%%%%%%%%%%%% 
ax2=subplot(2,3,4);

hold on
plot(xi_orig(1:5:end), f_orig(1:5:end), '-', 'LineWidth',2, 'color',c1,'DisplayName','SMP','MarkerSize',2);
plot(xi_orig1, f_orig1, lins, 'LineWidth',2, 'color',c2,'DisplayName','MP');
t_min = 200;
t_max = 1000;
idx = (xi_orig >= t_min) & (xi_orig <= t_max);
poly = polyfit(log10(xi_orig(idx)), log10(f_orig(idx)), 1);
y_fit = 10.^(polyval(poly, log10(xi_orig(idx))));
plot(xi_orig(idx), y_fit*0.1, '-.', 'LineWidth',1.5, 'color',c3);
text(24.246276965608793,1.66774147e-7,0,['slope:$\,$' num2str(round(poly(1),1))],'fontname','Timesnewroman','FontSize',13)
set(gca,'YScale','log','XScale','log')
% xlabel('FPT (s)')
% ylabel('PDF')
xlim([10 t])
ylim([10^-10 0.2])
% xticks([10 100])
yticks([1e-8  1e-5 1e-2])
xticks([ 1e2  1e4 ])
% xlabel('FPT (s)')
% ylabel('PDF')
set(gca, 'TickLength', [0.02, 0.025]);
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% xd=200
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
load("a6_FPT_20260609_200.mat")
%%%%
sur(1)=1-sum(isnan(FPT))/n;
sur1(1)=1-sum(isnan(FPT1))/n;
%%%%
ax3=subplot(2,3,2);

hold on;
%%%%
FPT_clean = FPT(~isnan(FPT));
FPT_clean1 = FPT1(~isnan(FPT1));
%%%%%%%%%% SMP
[f_log, xi_log] = ksdensity(log(FPT_clean), 'NumPoints', 1000, 'Support', 'positive','Bandwidth','plug-in');
% 
xi_orig = exp(xi_log);                       
f_orig = f_log ./ xi_orig;                   
%%%%%%%%%% MP
[f_log1, xi_log1] = ksdensity(log(FPT_clean1), 'NumPoints', 1000, 'Support', 'positive','Bandwidth','plug-in');
% 
xi_orig1 = exp(xi_log1);                     
f_orig1 = f_log1 ./ xi_orig1;                    
%
plot(xi_orig, f_orig, '-', 'LineWidth',2, 'color',c1,'DisplayName','SMP','MarkerSize',3);
plot(xi_orig1, f_orig1, lins, 'LineWidth',2, 'color',c2,'DisplayName','MP');
% legend('Box','off','Location','northeast')

% ylabel('PDF')
xlim([5 200])
ylim([0 0.035])
xticks([50 100 150])
% yticks(0:0.01:0.03)
set(gca,'YTickLabel',[])
title('$x_d=200\ \mu \rm m$','FontSize',14)
box on

%%%%%%%%%%%%% 
ax4=subplot(2,3,5);

hold on
plot(xi_orig(1:5:end), f_orig(1:5:end), '-', 'LineWidth',2, 'color',c1,'DisplayName','SMP','MarkerSize',2);
plot(xi_orig1, f_orig1, lins, 'LineWidth',2, 'color',c2,'DisplayName','MP');
t_min = 200;
t_max = 1000;
idx = (xi_orig >= t_min) & (xi_orig <= t_max);
poly = polyfit(log10(xi_orig(idx)), log10(f_orig(idx)), 1);
y_fit = 10.^(polyval(poly, log10(xi_orig(idx))));
plot(xi_orig(idx), y_fit*0.1, '-.', 'LineWidth',1.5, 'color',c3);
text(52.46937858857361,5.06295837e-7,0,['slope:$\,$' num2str(round(poly(1),1))],'fontname','Timesnewroman','FontSize',13)
set(gca,'YScale','log','XScale','log')

xlim([10 t])
ylim([10^-10 0.2])
yticks([1e-8  1e-5 1e-2])
xticks([ 1e2  1e4 ])
xlabel('First passage time (s)')
% ylabel('PDF')
set(gca,'YTickLabel',[])
set(gca, 'TickLength', [0.02, 0.025]);
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% xd=400
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
load("a6_FPT_20260609_400.mat")
%%%%
sur(1)=1-sum(isnan(FPT))/n;
sur1(1)=1-sum(isnan(FPT1))/n;
%%%%
ax5=subplot(2,3,3);

hold on;
%%%%
FPT_clean = FPT(~isnan(FPT));
FPT_clean1 = FPT1(~isnan(FPT1));
%%%%%%%%%% SMP
[f_log, xi_log] = ksdensity(log(FPT_clean), 'NumPoints', 1000, 'Support', 'positive','Bandwidth','plug-in');
% 
xi_orig = exp(xi_log);                     
f_orig = f_log ./ xi_orig;                   
%%%%%%%%%% MP
[f_log1, xi_log1] = ksdensity(log(FPT_clean1), 'NumPoints', 1000, 'Support', 'positive','Bandwidth','plug-in');
% 
xi_orig1 = exp(xi_log1);                      
f_orig1 = f_log1 ./ xi_orig1;                   
%
plot(xi_orig, f_orig, '-', 'LineWidth',2, 'color',c1,'DisplayName','SMP','MarkerSize',3);
plot(xi_orig1, f_orig1, lins, 'LineWidth',2, 'color',c2,'DisplayName','MP');
% legend('Box','off','Location','northeast')

% ylabel('PDF')
xlim([5 200])
ylim([0 0.035])
xticks([50 100 150])
% yticks(0:0.01:0.03)
set(gca,'YTickLabel',[])
title('$x_d=400\ \mu \rm m$','FontSize',14)
box on
ylabel('Linear','Units','normalized','Position',[1.17,0.52],'FontSize',13)
%%%%%%%%%%%%% 
ax6=subplot(2,3,6);

hold on
plot(xi_orig(1:5:end), f_orig(1:5:end), '-', 'LineWidth',2, 'color',c1,'DisplayName','SMP','MarkerSize',2);
plot(xi_orig1, f_orig1, lins, 'LineWidth',2, 'color',c2,'DisplayName','MP');
t_min = 200;
t_max = 1000;
idx = (xi_orig >= t_min) & (xi_orig <= t_max);
poly = polyfit(log10(xi_orig(idx)), log10(f_orig(idx)), 1);
y_fit = 10.^(polyval(poly, log10(xi_orig(idx))));
plot(xi_orig(idx), y_fit*0.1, '-.', 'LineWidth',1.5, 'color',c3);
text(22.84839851333908,0.000001119146522,0,['slope:$\,$' num2str(round(poly(1),1))],'fontname','Timesnewroman','FontSize',13)
set(gca,'YScale','log','XScale','log')

xlim([10 t])
ylim([10^-10 0.2])
yticks([1e-8  1e-5 1e-2])
xticks([ 1e2  1e4 ])
% xlabel('FPT (s)')
% ylabel('PDF')
set(gca,'YTickLabel',[])
set(gca, 'TickLength', [0.02, 0.025]);
ylabel('Log-log','Units','normalized','Position',[1.17,0.506],'FontSize',13)
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
text(-1.3,2.15,'(a)','Units','normalized','Interpreter','tex','FontName','Timesnewroman','FontSize',16)
text(-0.3,2.15,'(b)','Units','normalized','Interpreter','tex','FontName','Timesnewroman','FontSize',16)
text(0.7,2.15,'(c)','Units','normalized','Interpreter','tex','FontName','Timesnewroman','FontSize',16)
text(-1.3,0.88,'(d)','Units','normalized','Interpreter','tex','FontName','Timesnewroman','FontSize',16)
text(-0.3,0.88,'(e)','Units','normalized','Interpreter','tex','FontName','Timesnewroman','FontSize',16)
text(0.74,0.88,'(f)','Units','normalized','Interpreter','tex','FontName','Timesnewroman','FontSize',16)

set(ax1,'Position',[0.151,0.59,0.262,0.341])
set(ax2,'Position',[0.151,0.16,0.262,0.341])
set(ax3,'Position',[0.4124,0.59,0.262,0.341])
set(ax4,'Position',[0.4124,0.16,0.262,0.341])
set(ax5,'Position',[0.674,0.59,0.262,0.341])
set(ax6,'Position',[0.674,0.16,0.262,0.341])