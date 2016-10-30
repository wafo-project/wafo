%ITMKURS_LAB3 Script to computer exercises 3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Power Spectrum and Rainflow Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%%Power Spectrum and Rainflow Analysis
%%Power Spectrum
%

S=jonswap;
S=oscspec;
S=torsethaugen;
load deep.dat, S=dat2spec(deep); 

%Solution Problem1:

L=spec2mom(S,4);
f0=sqrt(L(2)/L(1))/2/pi
s=sqrt(L(1)) 
alfa=f0/(sqrt(L(3)/L(2))/2/pi)
wspecplot(S);
2*pi*f0

beta=3:0.1:5;
gam=5E-9;
dpl=gam*spec2dplus(S,beta);
Tfpl=(1./dpl)/3600/24  %in days
clf, plot(beta,Tfpl)

%
%%Simulation of $X$ and the damage intensity
%


T=1000/f0 % in seconds
gam=5E-9;
100*T*gam*spec2dplus(S,4.22) % in percent


max(S.w)/pi
dt=1./f0/60
S1=specinterp(S,dt);
[max(S1.w) pi/dt]
clf, plotspec(S), hold on, plotspec(S1,1,'r.')
clf, plot(beta,Tfpl)
hold on
for  i=1:5
   X=spec2sdat(S1,60000);
   tp=dat2tp(X);
   rfc=tp2rfc(tp);
   db=cc2dam(rfc,beta);
   plot(beta,(1000/f0)*(1./db)*(1/gam)/3600/24,'r.') % Tf in days
end

%
%%Theoretical computation of damage intensity
%


S=createspec;
w=levels([0 4 253]);
S.S=2*exp(-2*abs(w-2).^1.4);
S.S=S.S+9*exp(-8*(w-0.5).^2);
S.w=w; 
S.note=['S(w)=9*exp(-8(w-0.5)^2)+2*exp(-2|w-2|.^1.4)']
wspecplot(S)

L=spec2mom(S,4);
f0=sqrt(L(2)/L(1))/2/pi
dpl=gam*spec2dplus(S,beta);
Tfpl=(1./dpl)/3600/24  %in days
dt=1./f0/60
S1=specinterp(S,dt);
figure(1)
clf, plot(beta,Tfpl)
hold on
for  i=1:10
  X=spec2sdat(S1,60000);
  tp=dat2tp(X);
  rfc=tp2rfc(tp);
  db=cc2dam(rfc,beta);
  plot(beta,(1000/f0)*(1./db)*(1/gam)/3600/24,'r.') % Tf in days
end



a=sqrt(L(3)/L(2))/2/pi
s=sqrt(L(1))
help spec2cmat
paramu=[-4.5*s 4.5*s 46];
nit=2;
figure(2), clf
[frfc fMm]=spec2cmat(S,[],'rfc',[],paramu,nit);
hold on, plot(rfc(:,2),rfc(:,1),'.')
clf, pdfplot(fMm)
dg=a*cmat2dam(paramu,frfc.f,beta);
figure(1), plot(beta,(1./dg)*(1/gam)/3600/24,'g') % Tf in days



figure(2), clf
[frfc1 fMm1]=spec2cmat(S,[],'rfc',[],paramu,3);
dg1=a*cmat2dam(paramu,frfc1.f,beta);
figure(1)
plot(beta,(1./dg1)*(1/gam)/3600/24,'k') % Tf in days


help mctpsim
F{1,1}=fMm.f;
F{1,2}=fMm.f';
mctp=mctpsim(F,1000);
clf,plot(mctp(1:40)) 

frfc1=zeros(size(frfc.f));
frfc1=triu(frfc.f,5);
fMm1=iter(frfc1,fMm.f,20,0.001);
clf,contour(fMm.f,20) 
hold on,  contour(fMm1,20)
F{1,1}=fMm1;
F{1,2}=fMm1';
mctp1=mctpsim(F,1000);
clf, subplot(2,1,1)
plot(mctp(1:40))
subplot(2,1,2)
plot(mctp1(1:40))








