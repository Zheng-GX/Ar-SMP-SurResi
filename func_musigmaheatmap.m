%%%%%%%%%%
function out = func_musigmaheatmap(t,dt,n,iftrunc,varargin)
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
%%
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
Punif=@(n)2*pi*rand(n,1);
Dr=0.062; 
r0=15;
v0=20;
if iftrunc == 0
    %%%%%%%%
    lambda=exp(mu+sigma^2/2);
else
    %%%%%%%%
    lambda=func_meantrunc(mu,sigma,T);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hazard_f_exp = 1/lambda;
hazard_f = @(t) normpdf((log(t)-mu)/sigma) ./ (t * sigma .* (1 - normcdf((log(t)-mu)/sigma))); 
dimension=2;
%%%%%%%%%%%%%%%%%%%
x_grid = linspace(-5, 20, 500);
v_u_all=zeros(n,imax);
v_u1_all=zeros(n,imax);
for ih=1:2
    tauplus=zeros(n,1);
    tauminu=zeros(n,1);
    tautumble=zeros(n,1);
    b=ones(n,1);
    c=ones(n,1);
    d=ones(n,1);
    s=ones(n,1);
    s1b0=zeros(n,1);
    ia=ones(n,1);
    ib=ones(n,1);
    Ts=zeros(n,imax+1);
    Tb=zeros(n,imax+1);
    m=m_0*ones(n,1);
    a=a0*ones(n,1);
    Y=Y0*ones(n,1);
    theta=Punif(n); 
    x=zeros(n,imax+1);
    y=zeros(n,imax+1);
    for i=1:imax
        for j=1:n  
            if s(j)==1
                x(j,i+1)=x(j,i)+cos(theta(j,i))*v0*dt;
                y(j,i+1)=y(j,i)+sin(theta(j,i))*v0*dt;
                if s1b0(j)==1
                    theta(j,i+1) = wrapTo2Pi(theta(j,i)-v0/r0*dt + sqrt(2*Dr*dt) * randn); 
                else
                    theta(j,i+1) = wrapTo2Pi(theta(j,i) + sqrt(2*Dr*dt) * randn); 
                end
                %%%%%%%%%%%%%
                Epsilon=alpha*(m0-m(j,i))+log((1+L(x(j,i+1),y(j,i+1))/KI)/(1+L(x(j,i+1),y(j,i+1))/KA));
                a(j,i+1)=1/(1+exp(N*Epsilon));
                m(j,i+1)=m(j,i)+dt*(kR*(1-a(j,i))-kB*a(j,i));
                Y(j,i+1)=Ya*a(j,i);
                P1=0.2^-1*Y(j,i+1)^H/Y5^H;
                %累
                if L(x(j,i+1),y(j,i+1))>L(x(j,i),y(j,i))
                    tauplus(j,b(j))=tauplus(j,b(j))+dt;
                else
                    if L(x(j,i+1),y(j,i+1))<L(x(j,i),y(j,i))
                        tauminu(j,c(j))=tauminu(j,c(j))+dt;
                    end
                end
                if P1*dt>rand
                    s(j)=0;
                end
            else
                if s1b0(j)==1
                    theta(j,i+1)=wrapTo2Pi(theta(j,i)+TumbleAngleSampling1(1));% exp.surface
                else
                    theta(j,i+1)=wrapTo2Pi(theta(j,i)+(-1).^round(rand)*TumbleAngleSampling_Berg(1)); % exp.bulk
                end
                x(j,i+1)=x(j,i);
                y(j,i+1)=y(j,i);
                %%%%%%%%%%%%%
                Epsilon=alpha*(m0-m(j,i))+log((1+L(x(j,i+1),y(j,i+1))/KI)/(1+L(x(j,i+1),y(j,i+1))/KA));
                a(j,i+1)=1/(1+exp(N*Epsilon));
                m(j,i+1)=m(j,i)+dt*(kR*(1-a(j,i))-kB*a(j,i));
                Y(j,i+1)=Ya*a(j,i);
                P0=tau0^-1;
                tautumble(j,d(j))=tautumble(j,d(j))+dt;
                if P0*dt>rand
                    s(j)=1;
                    b(j)=b(j)+1;
                    c(j)=c(j)+1;
                    d(j)=d(j)+1;
                    tauplus(j,b(j))=0;
                    tauminu(j,c(j))=0;
                    tautumble(j,d(j))=0;
                end
            end
            if s1b0(j)==1
                Ts(j,ia(j))=Ts(j,ia(j))+dt;
                if ih==1
                    alpha_sb=hazard_f(Ts(j,ia(j)));
                else
                    alpha_sb=hazard_f_exp;
                end
                if rand<alpha_sb*dt
                    s1b0(j)=0;
                    ia(j)=ia(j)+1;
                end
            else
                Tb(j,ib(j))=Tb(j,ib(j))+dt;
                if rand<kb*dt
                    s1b0(j)=1;
                    ib(j)=ib(j)+1;
                end
            end
            tp(j)=sum(tauplus(j,:));
            tm(j)=sum(tauminu(j,:));
            tt(j)=sum(tautumble(j,:));
            if ih==1
                v_u_all(j,i)=v0*(tp(j)-tm(j))/(tp(j)+tm(j)+2*tt(j))/dimension;
            else
                v_u1_all(j,i)=v0*(tp(j)-tm(j))/(tp(j)+tm(j)+2*tt(j))/dimension;
            end
        end
    end
    mtp=mean(tp);
    mtm=mean(tm);
    mtt=mean(tt);
    mv=v0*(mtp-mtm)/(mtp+mtm+2*mtt)/dimension;
    if ih==1
        mvd=mv;
    else
        mvd1=mv; 
    end
end
eps = 1e-12;
% 
vd  = v_u_all(:, end);
vd1 = v_u1_all(:, end);
f_a  = ksdensity(vd,  x_grid, 'NumPoints', 500,'Bandwidth',0.2);
f1_a = ksdensity(vd1, x_grid, 'NumPoints', 500,'Bandwidth',0.2);
% 
f_a  = f_a + eps;
f1_a = f1_a + eps;
f  = f_a  / trapz(x_grid, f_a);
f1 = f1_a / trapz(x_grid, f1_a);

% 
integrand = f .* log(f ./ f1);
integrand(isnan(integrand) | isinf(integrand)) = 0;
KL_M2E = trapz(x_grid, integrand);

%
integrand = f1 .* log(f1 ./ f);
integrand(isnan(integrand) | isinf(integrand)) = 0;
KL_E2M= trapz(x_grid, integrand);
out.KL_M2E = KL_M2E;
out.KL_E2M = KL_E2M;
out.mvd = mvd;
out.mvd1 = mvd1;
out.vu=v_u_all;
out.vu1=v_u1_all;
end