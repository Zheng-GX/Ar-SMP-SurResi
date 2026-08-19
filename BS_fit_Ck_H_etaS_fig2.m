close all;clear;clc
%% data
surH    = [65, 80, 115, 160];     
surns10 = [0.471290828, 0.396886736, 0.345141876, 0.285099861];
%% parameter
global D t dt n mu sigma
t = 200;
dt = 0.01;             
imax = t / dt;         
n = 5000;     
D = 87;               
mu=2.36;
sigma=1.16;
%%
p0=12;  %kb~1/H or ~1/H^2
% 
diary('optim_log.txt');    
options = optimset('Display','iter','TolX',1e-6,'TolFun',1e-4);
pfit = fminsearch(@(p) sum((func_eta(p, surH) - surns10).^2), p0, options);

diary off;                 
 eta_pred = func_eta(pfit, surH);
 SSE = sum((eta_pred - surns10).^2);
 SST = sum((surns10 - mean(surns10)).^2);
 R2 = 1 - SSE / SST;
%% 
High = linspace(40,250,25);  
mps=func_eta(pfit, High);
timeStr = datestr(now, 'yyyymmdd_HHMMSS');
save(['a6_fit','_', timeStr '.mat']);
