close all;clear;clc
rng('shuffle');
%%
t=200;
dt=0.01;
n=5000;
High=linspace(30,250,20);
Gradient=[2 1 0.5];
%High=30;
%Gradient=0.1;
%%
for iH=1:length(High)
    H=High(iH);
    for iG=1:length(Gradient)
        G=Gradient(iG);
        resu=func_vd(t,dt,n,G,H);
        vdx_1(iH,iG)=resu.vd;
    end
end
% 
timeStr = datestr(now, 'yyyymmdd_HHMMSS');
save(['3G_H_vd_', timeStr, '.mat']);
