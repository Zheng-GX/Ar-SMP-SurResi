%%%%%%%%%%
function out = func_FPT(t,dt,n,mu,sigma,g,destination)
%%
N=6;
KI=18.2;
KA=3*10^3;
alpha=1.7;
m0=1; 
kR=0.005;
kB=kR;
a0=kR/(kR+kB);
Y5=3;
H=10;
Ya=Y5*(0.2/0.8)^(1/H)/a0;
Y0=Ya*a0;
tau0=0.2;
v0=20;
%
G=g*10^-3; 
L0=3*KI;
L=@(x,y)L0*exp(G*x);
Punif=@(n)2*pi*rand(n,1);
m_0=1-log((1+L0/KA)/(1+L0/KI))/alpha;
Dr=0.062;
r0=15;
p_bs=0.05;
lambda=exp(mu+sigma^2/2);
hazard_f_exp = 1/lambda;
hazard_f = @(t) normpdf((log(t)-mu)/sigma) ./ (t * sigma .* (1 - normcdf((log(t)-mu)/sigma))); 
FPT=NaN(n,1);  
FPT1=NaN(n,1);
%%
for ih=1:2
    x=zeros(n,1);
    y=zeros(n,1);
    s=ones(n,1);
    m=m_0*ones(n,1);
    a=a0*ones(n,1);
    Y=Y0*ones(n,1);
    theta=Punif(n);   
    s1b0=zeros(n,1);
    for j=1:n
        i=0;
        Ts=0;
        while x(j)<destination && i<t
            if s(j)==1
                x(j)=x(j)+cos(theta(j))*v0*dt;
                y(j)=y(j)+sin(theta(j))*v0*dt;
                if s1b0(j)==1
                    theta(j) = wrapTo2Pi(theta(j)-v0/r0*dt + sqrt(2*Dr*dt) * randn); 
                else
                    theta(j) = wrapTo2Pi(theta(j) + sqrt(2*Dr*dt) * randn); 
                end
                %%%%%%%%%%%%%
                Y(j)=Ya*a(j);
                Epsilon=alpha*(m0-m(j))+log((1+L(x(j),y(j))/KI)/(1+L(x(j),y(j))/KA));
                a(j)=1/(1+exp(N*Epsilon));
                m(j)=m(j)+dt*(kR*(1-a(j))-kB*a(j));
                P1=0.2^-1*Y(j)^H/Y5^H;
                if P1*dt>rand
                    s(j)=0;
                end
            else
                if s1b0(j)==1
                    theta(j)=wrapTo2Pi(theta(j)+TumbleAngleSampling1(1));
                else
                    theta(j)=wrapTo2Pi(theta(j)+(-1).^round(rand)*TumbleAngleSampling_Berg(1)); 
                end
                %%%%%%%%%%%%%
                Y(j)=Ya*a(j);
                Epsilon=alpha*(m0-m(j))+log((1+L(x(j),y(j))/KI)/(1+L(x(j),y(j))/KA));
                a(j)=1/(1+exp(N*Epsilon));
                m(j)=m(j)+dt*(kR*(1-a(j))-kB*a(j));
                P0=tau0^-1;
                if P0*dt>rand
                    s(j)=1;
                end
            end
            if s1b0(j)==1
                Ts=Ts+dt;
                if ih==1
                    alpha_sb=hazard_f(Ts);
                else
                    alpha_sb=hazard_f_exp;
                end
                if rand<alpha_sb*dt
                    s1b0(j)=0;
                    Ts=0;
                end
            else
                if rand<p_bs*dt
                    s1b0(j)=1;
                end
            end
            i=i+dt;
        end
        %%%%%%FPT
        if abs(i - t) > 1e-5
            if ih==1
                FPT(j)=i;
            else
                FPT1(j)=i;
            end
        end
    end
end
%%
mFPT = mean(FPT, 'omitnan');
mFPT1 = mean(FPT1, 'omitnan');
out.mFPT=mFPT;
out.mFPT1=mFPT1;
out.FPT=FPT;
out.FPT1=FPT1;
end