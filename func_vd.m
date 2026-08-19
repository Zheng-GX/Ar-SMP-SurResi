function out = func_vd(t,dt,n,g,High)
%%
imax=t/dt;
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
L0=3*KI;
Punif=@(n)2*pi*rand(n,1);
x=zeros(n,imax+1);
y=zeros(n,imax+1);
s=ones(n,1);
m_0=1-log((1+L0/KA)/(1+L0/KI))/alpha;
m=m_0*ones(n,1);
a=a0*ones(n,1);
Y=Y0*ones(n,1);
theta=Punif(n);  
Dr=0.062; 
mu = 2.36;          
sigma = 1.16;
dimension=2;
ia=ones(n,1);
ib=ones(n,1);
s1b0=zeros(n,1);
Ts=zeros(n,imax+1);
Tb=zeros(n,imax+1);
hazard_f = @(t) normpdf((log(t)-mu)/sigma) ./ (t * sigma .* (1 - normcdf((log(t)-mu)/sigma))); 
%%%%%%%%%%%%%%%%%%%%%%%
% kb=@(H)6*tau0*v0^2./H.^2;
C1=2.91;
kb=@(H)C1./H;
r0=15;
%%
p_bs=kb(High);
G=g*10^-3;
L=@(x,y)L0*exp(G*x);
tauplus=zeros(n,1);
tauminu=zeros(n,1);
tautumble=zeros(n,1);
b=ones(n,1);
c=ones(n,1);
d=ones(n,1);
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
            %
            if L(x(j,i+1),y(j,i+1))>L(x(j,i),y(j,i))
                tauplus(j,b(j))=tauplus(j,b(j))+dt;
            else
                if L(x(j,i+1),y(j,i+1))<L(x(j,i),y(j,i))
                    tauminu(j,c(j))=tauminu(j,c(j))+dt;
                end
            end
            % 
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
            alpha_sb=hazard_f(Ts(j,ia(j)));
            if rand<alpha_sb*dt
                s1b0(j)=0;
                ia(j)=ia(j)+1;
            end
        else
            Tb(j,ib(j))=Tb(j,ib(j))+dt;
            if rand<p_bs*dt
                s1b0(j)=1;
                ib(j)=ib(j)+1;
            end
        end
        tp(j)=sum(tauplus(j,:));
        tm(j)=sum(tauminu(j,:));
        tt(j)=sum(tautumble(j,:));
        % v_u(j,i)=v0*(tp(j)-tm(j))/(tp(j)+tm(j)+2*tt(j))/dimension;
    end
end
tp_m=mean(tp);
tm_m=mean(tm);
tt_m=mean(tt);
vdx_1=v0*(tp_m-tm_m)/(tp_m+tm_m+2*tt_m)/dimension;

out.vd=vdx_1;
end