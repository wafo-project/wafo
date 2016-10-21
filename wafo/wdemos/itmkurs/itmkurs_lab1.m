%ITMKURS_LAB1 Script to computer exercises 1
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%Analysis of Load Data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A Stochastic Load Process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%Measured data
%

load deep.dat
x = deep;
plot(x(:,1),x(:,2))
plot(x(1:1000,1),x(1:1000,2))


T=x(end,1)-x(1,1);

whos


tp = dat2tp(x);
plot(x(:,1),x(:,2),tp(:,1),tp(:,2),'.-')
axis([0 100 -20 20])

tp1 = dat2tp(x,1);
plot(x(:,1),x(:,2),tp(:,1),tp(:,2),tp1(:,1),tp1(:,2))
axis([0 100 -20 20])


lc = tp2lc(tp);
lc(:,2)=lc(:,2)/T;
plot(lc(:,1),lc(:,2))
semilogx(lc(:,2),lc(:,1))


m=mean(x(:,2));
f0 = interp1(lc(:,1),lc(:,2),m,'linear');
f0



extr0=length(tp)/2/T;
alfa=f0/extr0

%
%%Gaussian process as a model for the deep water data
%


plotnorm(x(:,2))



S = dat2spec(deep);
wspecplot(S);

S
plot(S.w,S.S)



lam = spec2mom(S,4); L0=lam(1); L2=lam(2); L4=lam(3);


f0=1/(2*pi)*sqrt(L2/L0)
ux = -20:0.1:20;
ricex = f0*exp(-ux.*ux./(2*L0));
plot(lc(:,1),lc(:,2),'-',ux,ricex,'--')
semilogx(lc(:,2),lc(:,1),'-',ricex,ux,'--')

  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cycle counts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------
%Cycle Counts
%---------------------

proc = x(1:500,:);
democc

RFC = tp2rfc(tp);
mM = tp2mm(tp);


subplot(1,2,1), ccplot(RFC)
subplot(1,2,2), ccplot(mM)


ampRFC = cc2amp(RFC);
ampmM = cc2amp(mM);
subplot(1,2,1), hist(ampRFC)
subplot(1,2,2), hist(ampmM)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SN-data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%Calculation of damage intensity}
%

beta=3.2; gam=5.5E-10;
d_beta=cc2dam(RFC,beta)/T;
time_fail=1/gam/d_beta/3600 %in hours of the specific storm


%
%%Additional exercises, Optional
%

load SN

plot(N,s,'o')
axis([0 14e5 5 35 ])
loglog(N,s,'o')


plotnorm(reshape(log(N),8,5))


[e0,beta0,s20] = snplot(s,N,12)
[e0,beta0,s20] = snplot(s,N,14)
e0    = 5.5361e-10
beta0 = 3.2286
s20   = 0.0604

%
%%Calculation of the 95\% quantile for the fatigue life time
%
%%Fatigue life distribution under variable random load
%

D0 = e0*cumsum((RFC(:,2)-RFC(:,1)).^beta0);
plot(D0)

      
beta = 3:0.1:8;
DRFC = cc2dam(RFC,beta);
dRFC = DRFC/T
plot(beta,dRFC)


help ftf
[t0,F0] = ftf(e0,cc2dam(RFC,beta0)/T,s20,0.5,1);

[t1,F1] = ftf(e0,cc2dam(RFC,beta0)/T,s20,0,1);
[t2,F2] = ftf(e0,cc2dam(RFC,beta0)/T,s20,5,1);
plot(t0,F0,t1,F1,t2,F2)


taRFC = exp(-1.96*sqrt(0.06))/e0./dRFC;
DmM = cc2dam(mM,beta);
dmM = DmM/T
tamM = exp(-1.96*sqrt(0.06))/e0./dmM;
plot(beta,taRFC,beta,tamM,'r')

%
%%Crack growth data}  %%%%%%%%%%%%%
%




clear
load virkler

plot(v(:,2),v(:,1))


plot(v(:,2:69),v(:,1),'b-')



N = alevel(v,15);

plot(N,ones(1,length(N)),'o')

a = nlevel(v,2e5);
plot(a,ones(1,length(a)),'o')


