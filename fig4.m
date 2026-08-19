%% 
close all; clear; clc;
set(0,'defaultAxesFontName', 'Times new roman','defaultAxesFontsize',17);
set(0,'defaultTextInterpreter','latex'); 
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultfigurecolor','w')
set(0,'defaultAxeslinewidth',1)

%% ================== 

ml = 0.09; mr = 0.075; % 
mb = 0.1; % 
cgap = 0.09; rgap = 0.18; % 
cbw = 0.015; % 
axcbg= 0.008; % 

wt=(1-ml-mr-2*cgap)/3;
w =wt-axcbg-cbw;

h1 = 0.32;        
h2 = 0.32;         
mt = 1 - (mb + h2 + rgap + h1);  
y1 = mb + h2 + rgap;
y2 = mb;
x0 = @(i) ml + (i-1)*(wt+cgap);



fig = figure('Position',[100.5,279,1349,592.5]);


%% %%%%%%%%%%%%%%%% 
%% (a) si=2
load('a6_LogxheatkbkRDKL_si2.mat')
ax = axes('Position',[x0(1), y1, w, h1]); hold on;
set(ax,'fontsize',17)

lambda=exp(0.1+3.^2/2);
rRs=1./LkR./lambda;
rbs=1./Lkb./lambda;
[X, Y] = meshgrid(rRs, rbs);
surface(X, Y, zeros(size(X)), KLend', ...
    'EdgeColor', 'none', 'FaceColor', 'texturemap');
% pcolor(X, Y, KLend');

axis tight;       
%
colormap(ax, othercolor('BuOr_8'))
cb = colorbar(ax);
cb.Position = [x0(1)+w+axcbg, y1, cbw, h1];
ax.Position = [x0(1), y1, w, h1];          
cb.Title.String = '$D_{\rm KL}$'; 
cb.Title.Interpreter = 'latex';   
cb.FontSize=17;
cb.Ticks=0.03:0.05:2;   


box on
ylabel('$\tau_b/\left \langle \tau _s \right \rangle$', 'Interpreter','latex', 'Rotation',90);
xlabel('$\tau_R/\left \langle \tau _s \right \rangle$', 'Interpreter','latex');
xlim([0.01 10])
ylim([0.125 0.55])
set(ax,'XScale','log')
xticks([0.01 0.1 1 10])
A=[0.5,0.15];
B=[10, 0.15];
Bx=5.3;
modif=0.02;
scatter([ A(1) B(1)],[ A(2) B(2)],'k','filled')
text(A(1)+modif,A(2)+modif,'\it A','FontSize',17,'FontWeight','bold')
text(Bx,B(2)+modif,'\it B','FontSize',17,'FontWeight','bold')
ax.Layer = 'top';
title(['$\sigma=' num2str(sigma) ,'$'])
tabc=[-0.08,1.176];
text(ax,tabc(1),tabc(2), '(a)', ...
         'Units', 'normalized', ...
         'VerticalAlignment', 'top', ...
         'FontSize', 18,'Interpreter','tex','FontName','Times new roman');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dcm = datacursormode(gcf);

% set(dcm, 'UpdateFcn', @customDataCursor);


%% (b) si=2.5
load('a6_LogxheatkbkRDKL_si2.5.mat')
ax = axes('Position',[x0(2), y1, w, h1]); hold on;
set(ax,'fontsize',17)

rRs=1./LkR./lambda;
rbs=1./Lkb./lambda;
[X, Y] = meshgrid(rRs, rbs);
surface(X, Y, zeros(size(X)), KLend', ...
    'EdgeColor', 'none', 'FaceColor', 'texturemap');
% pcolor(X, Y, KLend');

axis tight;       
colormap(ax, othercolor('BuOr_8'))
cb = colorbar(ax);
cb.Position = [x0(2)+w+axcbg, y1, cbw, h1];
ax.Position = [x0(2), y1, w, h1];
cb.Title.String = '$D_{\rm KL}$'; 
cb.Title.Interpreter = 'latex';   
cb.FontSize=17;
cb.Ticks=0.1:0.2:0.7;
% clim([0.5 inf])

box on
ylabel('$\tau_b/\left \langle \tau _s \right \rangle$', 'Interpreter','latex', 'Rotation',90);
xlabel('$\tau_R/\left \langle \tau _s \right \rangle$', 'Interpreter','latex');
xlim([0.01 10])
ylim([0.125 0.55])
set(ax,'XScale','log')
xticks([0.01 0.1 1 10])
scatter([ A(1) B(1)],[ A(2) B(2)],'k','filled')
text(A(1)+modif,A(2)+modif,'\it A','FontSize',17,'FontWeight','bold')
text(Bx,B(2)+modif,'\it B','FontSize',17,'FontWeight','bold')
ax.Layer = 'top';
title(['$\sigma=' num2str(sigma) ,'$'])
text(ax, tabc(1),tabc(2), '(b)', ...
         'Units', 'normalized', ...
         'VerticalAlignment', 'top', ...
         'FontSize', 18,'Interpreter','tex','FontName','Times new roman');

%% (c) si=3
load('a6_LogxheatkbkRDKL_si3.mat')
ax = axes('Position',[x0(3), y1, w, h1]); hold on;
set(ax,'fontsize',17)

rRs=1./LkR./lambda;
rbs=1./Lkb./lambda;
[X, Y] = meshgrid(rRs, rbs);
surface(X, Y, zeros(size(X)), KLend', ...
    'EdgeColor', 'none', 'FaceColor', 'texturemap');
% pcolor(X, Y, KLend');

axis tight;       
colormap(ax, othercolor('BuOr_8'))
cb = colorbar(ax);
cb.Position = [x0(3)+w+axcbg, y1, cbw, h1];
ax.Position = [x0(3), y1, w, h1];
cb.Title.String = '$D_{\rm KL}$'; 
cb.Title.Interpreter = 'latex';   
cb.FontSize=17;
cb.Ticks=0.2:0.4:10;

box on
ylabel('$\tau_b/\left \langle \tau _s \right \rangle$', 'Interpreter','latex', 'Rotation',90);
xlabel('$\tau_R/\left \langle \tau _s \right \rangle$', 'Interpreter','latex');
xlim([0.01 10])
ylim([0.125 0.55])
set(ax,'XScale','log')
xticks([0.01 0.1 1 10])
scatter([ A(1) B(1)],[ A(2) B(2)],'k','filled')
text(A(1)+modif,A(2)+modif,'\it A','FontSize',17,'FontWeight','bold')
text(Bx,B(2)+modif,'\it B','FontSize',17,'FontWeight','bold')
ax.Layer = 'top';
title(['$\sigma=' num2str(sigma) ,'$'])
text(ax, tabc(1),tabc(2), '(c)', ...
         'Units', 'normalized', ...
         'VerticalAlignment', 'top', ...
         'FontSize', 18,'Interpreter','tex','FontName','Times new roman');
%%%%%%%%%%%%% insert 
insetAx = axes('Position', [0.812453669384728,0.793248945147683,0.089065974796145,0.125105485232068]); 
Lsi=[2 2.2 2.5 2.8 3];

for iDmax=1:length(Lsi)
readname_a='KLmax_sigma';
readname_b=num2str(Lsi(iDmax));
load([readname_a readname_b '.mat'])
meanKLmax(iDmax)=KLmax_mean;
stdKLmax(iDmax)=KLmax_std;
semKLmax(iDmax)=KLmax_sem;
end
color1=[228,123,8]./255;
errorbar(Lsi, meanKLmax, stdKLmax, 'o', 'MarkerSize', 4, ...
    'MarkerFaceColor', color1,'MarkerEdgeColor',color1, 'Color', 'red')
xlim([1.8 3.2])
ylim([0 2.8])
yticks([1 2])
xlabel('$\sigma$')
ylabel('$D_{\rm KL}^{\rm max}$','Rotation',0)
ax = gca;

ax.XRuler.TickLabelGapOffset = -5; 

ax.YRuler.TickLabelGapOffset = -2; 

text(0.89, 1, '$\tau_b/\left \langle \tau _s \right \rangle$=0.15', 'Units', ...
    'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top','FontSize',12);
box off
%% %%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
blue1=[59,115,172]/255;yellow=[184,150,37]/255;
colors = [yellow;blue1]; 

%% (d) si=2
ax = axes('Position',[x0(1), y2, wt, h2]); hold on;
load('a6_A_vd_pdf_si2.mat');
lambda=round(exp(mu+sigma^2/2),1);
%%%%%%%%%%%v_u
fill([xi, fliplr(xi)], [f, zeros(size(f))], colors(1, :), ...
    'FaceAlpha', 0.6, 'EdgeColor', 'none', 'HandleVisibility', 'off');
% 
h_uA = plot(xi, f, 'Color', colors(1, :), 'LineWidth', 1.5, 'DisplayName','SMP, $A$');
%%%%%%%%%%%v_u_exp
%
h_u1A = plot(xi1, f1, '--', 'Color', colors(1, :),  'LineWidth', 1.5, 'DisplayName','MP, $A$');
% text(1.437307692307693,0.624625550660793,num2str(mvd,'%.2f'),'FontSize',17,'FontWeight','bold')
% text(0.641153846153846,0.878370044052864,num2str(mvd1,'%.2f'),'FontSize',17,'FontWeight','bold')
%%%%%%%%%%% v_u B
load('a6_B_vd_pdf_si2.mat');
fill([xi, fliplr(xi)], [f, zeros(size(f))], colors(2, :), ...
    'FaceAlpha', 0.6, 'EdgeColor', 'none', 'HandleVisibility', 'off');
h_uB = plot(xi, f, '-','Color', colors(2, :), 'LineWidth', 1.5, 'DisplayName','SMP, $B$');
h_u1B = plot(xi1, f1, '--','Color', colors(2, :),  'LineWidth', 1.5, 'DisplayName','MP, $B$');
xlabel('$v_d$ ($\mu$m/s)');
ylabel('PDF');
xlim([-0.7, 2.6]);
ylim([0 1.3])
% yticks([0 0.4 0.8])
box on;
legend([h_uA,h_u1A,h_uB,h_u1B],'Location', 'northeast', ...
    'FontSize', 13, 'Box', 'off');
title(['$\sigma=' num2str(sigma) ,'$'])
text(ax, tabc(1),tabc(2), '(d)', ...
         'Units', 'normalized', ...
         'VerticalAlignment', 'top', ...
         'FontSize', 18,'Interpreter','tex','FontName','Times new roman');

%% (e) si=2.5
ax = axes('Position',[x0(2), y2, wt, h2]); hold on;
load('a6_A_vd_pdf_si2.5.mat');
lambda=round(exp(mu+sigma^2/2),1);
%%%%%%%%%%%v_u
fill([xi, fliplr(xi)], [f, zeros(size(f))], colors(1, :), ...
    'FaceAlpha', 0.6, 'EdgeColor', 'none', 'HandleVisibility', 'off');
h_uA = plot(xi, f, 'Color', colors(1, :), 'LineWidth', 1.5, 'DisplayName','SMP, $A$');
h_u1A = plot(xi1, f1, '--', 'Color', colors(1, :),  'LineWidth', 1.5, 'DisplayName','MP, $A$');
%%%%%%%%%%% v_u B
load('a6_B_vd_pdf_si2.5.mat');
fill([xi, fliplr(xi)], [f, zeros(size(f))], colors(2, :), ...
    'FaceAlpha', 0.6, 'EdgeColor', 'none', 'HandleVisibility', 'off');
h_uB = plot(xi, f, '-','Color', colors(2, :), 'LineWidth', 1.5, 'DisplayName','SMP, $B$');
h_u1B = plot(xi1, f1, '--','Color', colors(2, :),  'LineWidth', 1.5, 'DisplayName','MP, $B$');
xlabel('$v_d$ ($\mu$m/s)');
ylabel('PDF');
xlim([-0.7, 2.6]);
ylim([0 1.3])
box on;
legend([h_uA,h_u1A,h_uB,h_u1B],'Location', 'northeast', ...
    'FontSize', 13, 'Box', 'off');
title(['$\sigma=' num2str(sigma) ,'$'])
text(ax, tabc(1),tabc(2), '(e)', ...
         'Units', 'normalized', ...
         'VerticalAlignment', 'top', ...
         'FontSize', 18,'Interpreter','tex','FontName','Times new roman');

%% (f) si=3
ax = axes('Position',[x0(3), y2, wt, h2]); hold on;
load('a6_A_vd_pdf_si3.mat');
lambda=round(exp(mu+sigma^2/2),1);
%%%%%%%%%%%v_u
fill([xi, fliplr(xi)], [f, zeros(size(f))], colors(1, :), ...
    'FaceAlpha', 0.6, 'EdgeColor', 'none', 'HandleVisibility', 'off');
h_uA = plot(xi, f, 'Color', colors(1, :), 'LineWidth', 1.5, 'DisplayName','SMP, $A$');
h_u1A = plot(xi1, f1, '--', 'Color', colors(1, :),  'LineWidth', 1.5, 'DisplayName','MP, $A$');
%%%%%%%%%%% v_u B
load('a6_B_vd_pdf_si3.mat');
fill([xi, fliplr(xi)], [f, zeros(size(f))], colors(2, :), ...
    'FaceAlpha', 0.6, 'EdgeColor', 'none', 'HandleVisibility', 'off');
h_uB = plot(xi, f, '-','Color', colors(2, :), 'LineWidth', 1.5, 'DisplayName','SMP, $B$');
h_u1B = plot(xi1, f1, '--','Color', colors(2, :),  'LineWidth', 1.5, 'DisplayName','MP, $B$');
xlabel('$v_d$ ($\mu$m/s)');
ylabel('PDF');
xlim([-0.7, 2.6]);
ylim([0 1.3])
box on;
legend([h_uA,h_u1A,h_uB,h_u1B],'Location', 'northeast', ...
    'FontSize', 13, 'Box', 'off');
title(['$\sigma=' num2str(sigma) ,'$'])
text(ax,tabc(1),tabc(2), '(f)', ...
         'Units', 'normalized', ...
         'VerticalAlignment', 'top', ...
         'FontSize', 18,'Interpreter','tex','FontName','Times new roman');
%% 
function txt = customDataCursor(~, event)
    
    h = event.Target;
   
    pos = event.Position;  
   
    xData = get(h, 'XData');
    yData = get(h, 'YData');
    cData = get(h, 'CData');

    
    [~, idxX] = min(abs(xData(1,:) - pos(1)));
    [~, idxY] = min(abs(yData(:,1) - pos(2)));
    colorVal = cData(idxY, idxX);

    
    txt = {sprintf('X: %.4f', pos(1)), ...
           sprintf('Y: %.4f', pos(2)), ...
           sprintf('D_{KL}: %.4f', colorVal)};
end