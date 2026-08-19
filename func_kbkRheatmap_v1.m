%%%%%%%%%%
function out = func_kbkRheatmap_v1(t,dt,n,iftrunc,varargin)
% 
p = inputParser;
% 
addParameter(p, 'kR', 0.005, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'kb', 0.05, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'mu', 0.1, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'sigma', 3, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'T', 400, @(x) isnumeric(x) && isscalar(x));
parse(p, varargin{:});

kR = p.Results.kR;
kb = p.Results.kb;
mu = p.Results.mu;
sigma = p.Results.sigma;
T  = p.Results.T;
% rng(12345);

imax=t/dt;
N=6;
KI=18.2;
KA=3*10^3;
alpha=1.7;
m0=1;
kB=kR;
a0=kR/(kR+kB);
Y5=3;
H=10;
tau0=0.2;
Ya=Y5*(tau0/0.8)^(1/H)/a0;
Y0=Ya*a0;
%%%%%%%%%%%%%%%%%
L0=3*KI;
g=2;
G=g*10^-3; 
L=@(x,y)L0*exp(G*x);
%%%%%%%%%%%%%%%%%%%
m_0=1-log((1+L0/KA)/(1+L0/KI))/alpha;
Dr=0.062;
r0=15;
v0=20;
if iftrunc == 0
    %%%%%%%%
    lambda=exp(mu+sigma^2/2);
else
    %%%%%%%%
    func_meantrunc =@(mu,sigma,T) exp(mu + sigma.^2/2) .* ...
        (normcdf((log(T) - mu - sigma.^2) ./ sigma) ./ normcdf((log(T) - mu) ./ sigma));
    lambda=func_meantrunc(mu,sigma,T);
end
hazard_f_exp = 1/lambda;
hazard_f = @(t) normpdf((log(t)-mu)/sigma) ./ (t * sigma .* (1 - normcdf((log(t)-mu)/sigma))); 
dimension=2;
%%%%%%%%%%%%%%%%%%%

v_final_M = zeros(n,imax);   % MP
v_final_E = zeros(n,imax);   % exp
v_final   = zeros(n,imax);
for ih=1:2
    % 
    x = zeros(n,2);   y = zeros(n,2);
    theta = zeros(n,2);
    a = a0 * ones(n,2);
    m = m_0 * ones(n,2);
    Y = Y0 * ones(n,2);
    s = ones(n,2);        
    s1b0 = zeros(n,2);   
    %
    theta(:,1) = 2*pi*rand(n,1);

    % 
    tp_total = zeros(n,1);   tm_total = zeros(n,1);   tt_total = zeros(n,1);
    curr_Ts = zeros(n,1);    curr_Tb = zeros(n,1);
    
    for i=1:imax
        for j=1:n  
            % 
            x_old = x(j,1);   y_old = y(j,1);
            theta_old = theta(j,1);
            a_old = a(j,1);   m_old = m(j,1);   
            s_old = s(j,1);   s1b0_old = s1b0(j,1);

            if s_old ==1
                x_new = x_old + cos(theta_old)*v0*dt;
                y_new = y_old + sin(theta_old)*v0*dt;
                L_old = L(x_old, y_old);
                L_new = L(x_new, y_new);
                if s1b0_old == 1
                    theta_new = wrapTo2Pi(theta_old - v0/r0*dt + sqrt(2*Dr*dt)*randn);
                else
                    theta_new = wrapTo2Pi(theta_old + sqrt(2*Dr*dt)*randn);
                end
                %%%%%%%%%%%%%
                Epsilon = alpha*(m0 - m_old) + log((1+L_new/KI)/(1+L_new/KA));
                a_new = 1/(1+exp(N*Epsilon));
                m_new = m_old + dt*(kR*(1-a_old) - kB*a_old);   
                Y_new = Ya * a_new;
                %
                if L_new > L_old
                    tp_total(j) = tp_total(j) + dt;
                elseif L_new < L_old
                    tm_total(j) = tm_total(j) + dt;
                end
                % 
                P1 = 0.2^-1 * Y_new^H / Y5^H;
                if P1*dt>rand
                    s_new = 0;   % tumble
                else
                    s_new=1;
                end
            else
                % 
                x_new = x_old;   y_new = y_old;
                % 
                if s1b0_old == 1
                    theta_new = wrapTo2Pi(theta_old + TumbleAngleSampling1(1));
                else
                    theta_new = wrapTo2Pi(theta_old + (-1)^round(rand)*TumbleAngleSampling_Berg(1));
                end
                %%%%%%%%%%%%%
                L_cur = L(x_new, y_new);
                Epsilon = alpha*(m0 - m_old) + log((1+L_cur/KI)/(1+L_cur/KA));
                a_new = 1/(1+exp(N*Epsilon));
                m_new = m_old + dt*(kR*(1-a_old) - kB*a_old);
                Y_new = Ya * a_new;
                
                tt_total(j) = tt_total(j) + dt;
                P0 = 1/tau0;
                if P0*dt>rand
                    s_new = 1;   % run
                else
                    s_new = 0;
                end
            end
            if s1b0_old == 1
                curr_Ts(j) = curr_Ts(j) + dt;
                if ih==1
                    alpha_sb = hazard_f(curr_Ts(j));
                else
                    alpha_sb=hazard_f_exp;
                end
                if rand<alpha_sb*dt
                    s1b0_new = 0;
                    curr_Ts(j) = 0;
                else
                    s1b0_new=1;
                end
            else
                curr_Tb(j) = curr_Tb(j) + dt;
                if rand<kb*dt
                    s1b0_new = 1;
                    curr_Tb(j) = 0;
                else
                    s1b0_new=0;
                end
            end
            x(j,2) = x_new;   y(j,2) = y_new;   theta(j,2) = theta_new;
            a(j,2) = a_new;   m(j,2) = m_new;   Y(j,2) = Y_new;
            s(j,2) = s_new;   s1b0(j,2) = s1b0_new;
            
        end
        x(:,1) = x(:,2);   y(:,1) = y(:,2);   theta(:,1) = theta(:,2);
        a(:,1) = a(:,2);   m(:,1) = m(:,2);   Y(:,1) = Y(:,2);
        s(:,1) = s(:,2);   s1b0(:,1) = s1b0(:,2);

        v_final(:,i) = v0 * (tp_total - tm_total) ./ (tp_total + tm_total + 2*tt_total) / dimension;
    end

    % 
    mtp = mean(tp_total);
    mtm = mean(tm_total);
    mtt = mean(tt_total);
    mv = v0 * (mtp - mtm) / (mtp + mtm + 2*mtt) / dimension;
    % 
    if ih == 1
        mvd = mv;
        v_final_M = v_final;
    else
        mvd1 = mv;
        v_final_E = v_final;
    end
end
%% ========== 
tKL=1/4*t;
Ntail = round(tKL/dt);
idx0 = imax - Ntail + 1;

v_tail_M = v_final_M(:, idx0:imax); 
v_tail_E = v_final_E(:, idx0:imax);

KL_t = zeros(Ntail,1);
KL_t_rev = zeros(Ntail,1);

x_grid = linspace(-5,20,500);
eps0 = 1e-12;

for k = 1:Ntail
    vd  = v_tail_M(:,k);
    vd1 = v_tail_E(:,k);

    fa  = ksdensity(vd,  x_grid, 'NumPoints', 500, 'Bandwidth', 0.2);
    f1a = ksdensity(vd1, x_grid, 'NumPoints', 500, 'Bandwidth', 0.2);

    fa  = fa  + eps0;
    f1a = f1a + eps0;

    f  = fa  / trapz(x_grid, fa);
    f1 = f1a / trapz(x_grid, f1a);

    integrand = f .* log(f ./ f1);
    integrand(isnan(integrand) | isinf(integrand)) = 0;
    KL_t(k) = trapz(x_grid, integrand);

    integrand = f1 .* log(f1 ./ f);
    integrand(isnan(integrand) | isinf(integrand)) = 0;
    KL_t_rev(k) = trapz(x_grid, integrand);
end

KL_M2E = mean(KL_t);
KL_E2M = mean(KL_t_rev);
KL_M2E_end = KL_t(end);
KL_E2M_end = KL_t_rev(end);

out.KL_M2E = KL_M2E;
out.KL_E2M = KL_E2M;
out.KL_M2E_end = KL_M2E_end;
out.KL_E2M_end = KL_E2M_end;

out.KL_t = KL_t;
out.KL_t_rev = KL_t_rev;
out.mvd = mvd;
out.mvd1 = mvd1;
end