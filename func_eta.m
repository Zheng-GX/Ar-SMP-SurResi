%%
function mps = func_eta(pfit, High)
rng(0)
global D t dt n mu sigma
% mu=pfit(1);
% sigma=pfit(2);
Ck=pfit(1);
imax=t/dt;
% kb=@(H)Ck*D./H.^2;
kb=@(H)Ck./H;
hazard_f = @(t) normpdf((log(t)-mu)/sigma) ./ (t * sigma .* (1 - normcdf((log(t)-mu)/sigma))); 
for ipbs=1:length(High)
    p_bs(ipbs)=kb(High(ipbs));
    s1b0=zeros(n,1);
    sur1=s1b0;
    ia=ones(n,1);
    Ts=zeros(n,imax+1);
    for i=1:imax
        for j=1:n   
            if s1b0(j)==1
                Ts(j,ia(j))=Ts(j,ia(j))+dt;
                alpha_sb=hazard_f(Ts(j,ia(j)));
                if rand<alpha_sb*dt
                    s1b0(j)=0;
                    ia(j)=ia(j)+1;
                end
            else
                if rand<p_bs(ipbs)*dt
                    s1b0(j)=1;
                end
            end
        end
        sur1(:,i+1)=s1b0;
    end
    ss=sum(sur1,1);
    pstb=ss/n;
    mps(ipbs)=mean(pstb(100/dt:t/dt));
end
end