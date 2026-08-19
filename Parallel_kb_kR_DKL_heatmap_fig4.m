%% BS_kb_kR_DKL_parallel.m
close all; clear; clc;

% =========
numCores = str2double(getenv('SLURM_CPUS_PER_TASK'));
if isnan(numCores) || numCores < 1
    numCores = feature('numcores');
    fprintf('未检测到 SLURM_CPUS_PER_TASK，使用本地核心数：%d\n', numCores);
else
    fprintf('SLURM 分配的核心数：%d\n', numCores);
end

% ========= 
MAX_WORKERS = 48;  
if numCores > MAX_WORKERS
    fprintf('核心数 %d 超过上限 %d，将使用 %d\n', numCores, MAX_WORKERS, MAX_WORKERS);
    numCores = MAX_WORKERS;
end
fprintf('计划启动 %d 个 worker\n', numCores);

% =========
c = parcluster('local');
fprintf('修改前 c.NumWorkers = %d\n', c.NumWorkers);

if numCores > c.NumWorkers
    fprintf('需要 %d 个 worker，将临时设置 c.NumWorkers = %d\n', numCores, numCores);
    c.NumWorkers = numCores;   
end

% =========
try
    parpool(c, numCores); 
    pool = gcp('nocreate');
    fprintf('并行池已启动，实际 worker 数量：%d\n', pool.NumWorkers);
    useParallel = true;
catch ME
    warning('并行池启动失败，切换到串行模式。错误：%s', ME.message);
    useParallel = false;
end

%% ========= 
rng('shuffle');

t  = 800;
dt = 0.01;
n  = 5000;        
iftrunc=0;
T=t;
lambda=exp(0.1+3.^2/2);
%lambda=50

sigma=2.8
mu=log(lambda)-sigma.^2/2;


rbs=linspace(0.12, 0.6, 20); 
rRs = logspace(log10(0.01), log10(10), 20);

%rRs=linspace(0.07, 2.5, 25); 


Lkb =1./(rbs.*lambda);
LkR =1./(rRs.*lambda);

% Lkb = [linspace(0.0005,0.025,8) ...
%   linspace(0.03,0.08,8)];
% LkR = [linspace(0.0005,0.015,4) ...
%    linspace(0.02,0.1,12) ...
%    linspace(0.11,0.2,6)];

NkR = length(LkR);
Nkb = length(Lkb);

%% ========= 
KL_M2E     = zeros(NkR, Nkb);
KLend      = zeros(NkR, Nkb);
mvd        = zeros(NkR, Nkb);
mvd1       = zeros(NkR, Nkb);
M_m_SMP    = zeros(NkR, Nkb);
re_M_m_SMP = zeros(NkR, Nkb);

%% ========= 
parfor ikR = 1:NkR
    kR = LkR(ikR);

    % 
    KL_row     = zeros(1, Nkb);
    Kend_row   = zeros(1, Nkb);
    mvd_row    = zeros(1, Nkb);
    mvd1_row   = zeros(1, Nkb);
    Mm_row     = zeros(1, Nkb);
    reMm_row   = zeros(1, Nkb);

    for ikb = 1:Nkb
        kb = Lkb(ikb);


        resu = func_kbkRheatmap_v1(t, dt, n, iftrunc, 'kR', kR, 'kb', kb,'mu',mu,'sigma',sigma,'T',T); % iftrunc=1表示有限窗口, 0表示理论平均值

        
        KL_row(ikb)   = resu.KL_M2E;
        Kend_row(ikb) = resu.KL_M2E_end;
        mvd_row(ikb)  = resu.mvd;
        mvd1_row(ikb) = resu.mvd1;
        Mm_row(ikb)   = resu.mvd - resu.mvd1;
        reMm_row(ikb) = (resu.mvd - resu.mvd1) ./ resu.mvd;
    end

    
    KL_M2E(ikR,:)     = KL_row;
    KLend(ikR,:)      = Kend_row;
    mvd(ikR,:)        = mvd_row;
    mvd1(ikR,:)       = mvd1_row;
    M_m_SMP(ikR,:)    = Mm_row;
    re_M_m_SMP(ikR,:) = reMm_row;
end

%% 
timeStr = datestr(now, 'yyyymmdd_HHMMSS');
save(['a6_LogxheatkbkRDKL_','si',num2str(sigma),'_', timeStr, '.mat'], ...
     't','dt','n','Lkb','LkR', ...
     'KL_M2E','mvd','mvd1','M_m_SMP','re_M_m_SMP','mu','sigma','iftrunc','T','KLend');


delete(gcp('nocreate'));