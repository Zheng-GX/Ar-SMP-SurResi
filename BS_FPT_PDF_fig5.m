close all;clear;clc
%%
t=100000;
dt=0.01;
n=5000;
mu=0.1;
sigma=3;
g=2;
destination=400;
resu = func_FPT(t,dt,n,mu,sigma,g,destination);
FPT=resu.FPT;
FPT1=resu.FPT1;
%%
 timeStr = datestr(now, 'yyyymmdd_HHMMSS');
 save(['a6_FPT_',timeStr,'.mat'],'FPT','FPT1','dt','n','mu','sigma','t','g','destination')
