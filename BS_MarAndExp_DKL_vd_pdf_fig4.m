close all;clear;clc
%%
t=800;
dt=0.01;
imax=t/dt;
tspan=0:dt:t;
%
n=5000;
iftrunc=0;
T=t;
%%
lambda=exp(0.1+3.^2/2);
sigma=3
mu=log(lambda)-sigma.^2/2;

%A
%kb= 0.067;
%kR= 0.02;

%B
kb= 0.067;
kR= 0.001;


resu=func_musigmaheatmap(t,dt,n,iftrunc,'mu',mu,'sigma',sigma,'T',T,'kR',kR,'kb',kb); 
v_u=resu.vu;
v_u1=resu.vu1;
mvd = resu.mvd;
mvd1 = resu.mvd1;

[f, xi] = ksdensity(v_u(:, end), 'NumPoints', 500,'Bandwidth',0.2);
[f1, xi1] = ksdensity(v_u1(:, end), 'NumPoints', 500,'Bandwidth',0.2);
%%
timeStr = datestr(now, 'yyyymmdd_HHMMSS');
save(['a6_B_vd_pdf_','si',num2str(sigma),'_',  timeStr, '.mat'],'t','dt','n','T','iftrunc','xi','f','xi1','f1','mu','sigma','kb','kR','mvd','mvd1');
