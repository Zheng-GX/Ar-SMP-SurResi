close all;clear;clc
%
set(0,'defaultAxesFontName', 'Times new roman','defaultAxesFontsize',22);
set(0,'defaultTextInterpreter','latex','defaultLegendInterpreter','latex'); 
set(0,'defaultfigurecolor','w','defaultAxeslinewidth',1.5)
set(0,'defaultfigureposition',[500,300,450,320]);
blue=[0.1216, 0.4667, 0.7059];red=[0.8392, 0.1529, 0.1569];
%% H-etaS
figure('Position',[776,365,434,394])
set(gca,'Position',[0.21,0.23,0.67,0.74])
hold on

%
% load("a6_fit_12D_H2_20260408_111060.mat")
% plot(High,mps,'--s','Color','k','MarkerFaceColor','k','MarkerSize',7,'displayname','$k_b=12D/H^2$')

%
load("a6_fit_1_H_20260409_094821.mat")
plot(High,mps,'-o','Color','k','MarkerFaceColor','k','MarkerSize',7,'displayname','$k_b=2.91/H$')
xlabel('$H\,(\mu\rm m)$')
ylabel('$\eta_{\rm{s}}$','Rotation',0)

%
plot(surH,surns10,'^','Color',red,'MarkerFaceColor',red,'MarkerSize',12,'displayname','Exp. (Wei 2025)')
box on
legend('FontSize',18,'LineWidth',1,'Box','off','Location','northeast')
xlim([30 250])
ylim([0.18 0.6])
yticks([ 0.3  0.5])
xticks([50 150 250])
%% H-vd
figure('Position',[776,365,434,394])
set(gca,'Position',[0.21,0.23,0.67,0.74])
hold on
%
load("a6_Hvd_1_H_20260409_101037.mat")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(High,vdx,'-o','Color','k','MarkerFaceColor','k','MarkerSize',7,'displayname','$k_b=2.91/H$')
xlabel('$H\,(\mu\rm m)$')
ylabel('$\left \langle v_{d} \right \rangle (\mu \rm m/s)$')

%
plot(70,1.2,'^','Color',red,'MarkerFaceColor',red,'MarkerSize',12,'DisplayName','Exp. (Grognot 2021)')
box on
legend('FontSize',18,'LineWidth',1,'Box','off','Location','northwest')


xlim([30 250])
ylim([0.8 3.1])
xticks([50 150 250])
yticks([ 1 2 3])   