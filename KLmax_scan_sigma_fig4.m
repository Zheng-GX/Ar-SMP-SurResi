
close all; clear; clc;

%% ====================== 
sigma = 2;          
% ----------------------------------------------------------------------
% 其余参数全部固定：
sigma_ref = 3;                            
mu_fix    = 0.1;                         
lambda = exp(mu_fix + sigma_ref^2/2);     
mu     = log(lambda) - sigma^2/2;       

iftrunc = 0;          
T  = 800;           

rbs  = 0.15;                       
rRs  = linspace(0.1, 1, 10);        
nRep = 20;                          

% 
t  = 800;      
dt = 0.01;     
n  = 5000;     
% ======================================================================

% 
Lkb = 1/(rbs*lambda);              
LkR = 1./(rRs*lambda);            
NkR = numel(LkR);

fprintf('sigma=%.3f, mu=%.4f, lambda=%.3f\n', sigma, mu, lambda);
fprintf('rbs=%.3f  ->  kb=%.5f\n', rbs, Lkb);
fprintf('rRs = %s\n', mat2str(rRs, 4));
fprintf('kR  = %s\n', mat2str(LkR, 5));

%% ======================
numCores = str2double(getenv('SLURM_CPUS_PER_TASK'));
if isnan(numCores) || numCores < 1
    numCores = feature('numcores');
end
MAX_WORKERS = 48;
if numCores > MAX_WORKERS
    numCores = MAX_WORKERS;
end
c = parcluster('local');
if numCores > c.NumWorkers
    c.NumWorkers = numCores;   
end
try
    parpool(c, numCores);
    fprintf('并行池已启动：%d workers\n', gcp('nocreate').NumWorkers);
catch ME
    warning('并行池启动失败，串行运行：%s', ME.message);
end

%% ============ 
rng('shuffle');
seeds = randi([1, 2^31-1], nRep, 1);  

KLmax_all  = zeros(nRep,1);   
idxmax_all = zeros(nRep,1);   
KLend_all  = zeros(nRep,NkR); 
KLmean_all = zeros(nRep,NkR); 

parfor irep = 1:nRep
    rng(seeds(irep), 'twister');
    Kend_row  = zeros(1,NkR);
    Kmean_row = zeros(1,NkR);
    for ikR = 1:NkR
        kR = LkR(ikR);
        resu = func_kbkRheatmap_v1(t, dt, n, iftrunc, ...
                   'kR', kR, 'kb', Lkb, 'mu', mu, 'sigma', sigma, 'T', T);
        Kend_row(ikR)  = resu.KL_M2E_end;  
        Kmean_row(ikR) = resu.KL_M2E;    
    end
    [mx, im] = max(Kend_row);             
    KLmax_all(irep)  = mx;
    idxmax_all(irep) = im;
    KLend_all(irep,:)  = Kend_row;
    KLmean_all(irep,:) = Kmean_row;
end

%% ====================== 
KLmax_mean = mean(KLmax_all);          
KLmax_std  = std(KLmax_all);          
KLmax_sem  = KLmax_std/sqrt(nRep);    

idx_mode     = mode(idxmax_all);      
rRs_at_mode  = rRs(idx_mode);         
rRs_max_mean = mean(rRs(idxmax_all)); 

fprintf('\n===== 结果汇总 (sigma=%.2f) =====\n', sigma);
fprintf('单次 KLmax（%d 次）：%s\n', nRep, mat2str(KLmax_all', 4));
fprintf('KLmax 位置 idx（%d 次）：%s\n', nRep, mat2str(idxmax_all', 4));
fprintf('KLmax 均值        = %.4f\n', KLmax_mean);
fprintf('KLmax 标准差 std  = %.4f   <- 误差棒（可选）\n', KLmax_std);
fprintf('KLmax 标准误 sem  = %.4f   <- 误差棒（推荐）\n', KLmax_sem);
fprintf('峰值最频繁位置 idx= %d (rRs=%.3f)\n', idx_mode, rRs_at_mode);

%% ====================== 
timeStr = datestr(now, 'yyyymmdd_HHMMSS');
matFile = ['KLmax_sigma', num2str(sigma), '_', timeStr, '.mat'];
save(matFile, 'sigma','mu','lambda','rbs','rRs','Lkb','LkR', ...
     'nRep','t','dt','n','iftrunc','T', ...
     'KLmax_all','idxmax_all','KLend_all','KLmean_all', ...
     'KLmax_mean','KLmax_std','KLmax_sem','idx_mode','rRs_at_mode','rRs_max_mean');
fprintf('已保存：%s\n', matFile);

csvFile = 'KLmax_vs_sigma.csv';
if ~exist(csvFile, 'file')
    fid = fopen(csvFile, 'w');
    fprintf(fid, 'sigma,KLmax_mean,KLmax_std,KLmax_sem,idx_mode,rRs_at_mode\n');
else
    fid = fopen(csvFile, 'a');
end
fprintf(fid, '%.4f,%.6f,%.6f,%.6f,%d,%.4f\n', ...
        sigma, KLmax_mean, KLmax_std, KLmax_sem, idx_mode, rRs_at_mode);
fclose(fid);
fprintf('已追加：%s\n', csvFile);

delete(gcp('nocreate'));
