close all;clear;clc
%
set(0,'defaultAxesFontName', 'Times new roman','defaultAxesFontsize',22);
set(0,'defaultTextInterpreter','latex','defaultLegendInterpreter','latex'); 
set(0,'defaultfigurecolor','w','defaultAxeslinewidth',1.5)
set(0,'defaultfigureposition',[500,300,450,320]);
%% H-vd
close all

% 
files = dir('3G_H_vd_*.mat');
Nfiles = length(files);


% 
vdx_cell = cell(1, Nfiles);

for i = 1:Nfiles
    load(files(i).name);
    vdx_cell{i} = vdx_1;       
end

% 
vdx_3d = cat(3, vdx_cell{:});

% 
mean_vdx = mean(vdx_3d, 3);
std_vdx = std(vdx_3d, 0, 3); 
% 

colors=[124,77,121
    91,181,172
    216,179,101]/255;
%
figure('Position',[776,365,434,394])
set(gca,'Position',[0.21,0.23,0.67,0.74])
hold on
MS=8;
ls={'o','s','^'};

for k = 1:length(Gradient)
p(k)=plot(High,mean_vdx(:,k),['-' ls{k}],'MarkerSize',MS,'Color',colors(k,:),'MarkerFaceColor',colors(k,:),...
    'displayname',['$G = \mathrm{', num2str(Gradient(k)), '\,mm}^{-1}$']);

% errorbar(High, mean_vdx(:,k), std_vdx(:,k), ...
%         'Color', colors(k,:), 'LineStyle', 'none', 'CapSize', 6, ...
%         'Marker', 'none', 'HandleVisibility', 'off');  % 不显示在 legend
end

xlabel('$H\,(\mu \rm m)$')
ylabel('$\left \langle v_{d} \right \rangle (\mu \rm m/s)$')
box on
legend('fontsize',17,'box','off','LineWidth',1, 'Interpreter', 'latex', 'FontName', 'Times New Roman')
% ylim([0 3])
ylim([0 2.5])
xlim([30 250])
yticks([0 1 2])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% 
files1 = dir('3H_kR_vd_*.mat');
Nfiles1 = length(files1);
if Nfiles1 == 0
    error('未找到任何 3H_kR_vd_*.mat 文件');
end
% 
vdx_cell1 = cell(1, Nfiles1);
for i = 1:Nfiles1
    load(files1(i).name);
    vdx_cell1{i} = vdx_1;       
end
% 
vdx_3d1 = cat(3, vdx_cell1{:});

% 
mean_vdx1 = mean(vdx_3d1, 3);
std_vdx1 = std(vdx_3d1, 0, 3); 

c1= [124,77,121
    91,181,172
    216,179,101
    41,109,119
    ]/255;
ls1={'v','d','h','o','s'};
figure('Position',[776,365,434,394])
set(gca,'Position',[0.21,0.23,0.67,0.74])
hold on
for i=1:length(High)
plot(kR_range,mean_vdx1(i,:),['-' ls1{i}],'MarkerSize',MS,'Color',c1(i,:),'MarkerFaceColor',c1(i,:),...
    'displayname',['$H=',num2str(High(i)),'\,\mu\rm m$'])

% errorbar(kR_range, mean_vdx1(i,:), std_vdx1(i,:), ...
%         'Color',c1(i,:), 'LineStyle', 'none', 'CapSize', 6, ...
%         'Marker', 'none', 'HandleVisibility', 'off');  % 不显示在 legend
end
xlabel('$k_{R}\,\rm (s^{-1})$')
ylabel('$\left \langle v_{d} \right \rangle (\mu \rm m/s)$')
box on
legend('fontsize',17,'box','off','LineWidth',1)
% ylim([1.52 4.5])
ylim([0 3.2])
xlim([0.0011 0.5])
% xlim([0.001 0.5])
xticks([ 0.01 0.1])
set(gca,'XScale','log') 