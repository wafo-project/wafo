function [mu,mu2,mu1,muG,beta,cc,fu,Fu]= chitwo2lc_sp(gamma,beta,S12,S22,ds0)
%CHITWO2LC_SP Saddlepoint approximation of crossing intensity, noncentral Chi^2 process 
%
% The noncentral Chi^2 process is defined in IR-(27).
%           X(t)=sum beta_j Z_j(t) + gamma_j Z_j(t)^2 
%
%  CALL:  [mu,fu,Fu]= chitwo2lc_sp(gamma,beta,S12,S22,ds);
%  
%       mu   = three column vector containing: levels in the first column, 
%              the saddlepoint approximation of crossing intensitity for quadratic sea  
%              in the second column and for Gaussian sea in the third column.
%              Note that the levels $u$ are not equally spaced. If the spacing 
%              in the grid is not fine enough, choose smaller value of parameter ds0. 
%       fu   = the pdf structure: f.x contains levels u, f.f is a two column matrix containing
%              the Danniel's-saddlepoint approximation of the density (pdf) of the quadratic 
%              sea in the first column and pdf of the Gaussian sea in the second.  
% gamma,beta = the vector containing constants defining the process
%              approximation of the cdf of the quadratic sea.  
%     S12    = covariance between vector Z(t) and Z'(t), where Z(t)=(Z_1(t),...,Z_n(t)).
%     S22    = covariance between vector Z'(t) and Z'(t)., where Z(t)=(Z_1(t),...,Z_n(t)).
%     ds0    = parameter defining the levels: default value is ds0=0.1.
%
% Example:
%   S=jonswap;
%   [gamma beta S12 S22]=dirsp2chitwo(S.S,S.w);
%   [mu, fu, Fu]=chitwo2lc_sp(gamma,beta,S12,S22); 
%   semilogy(mu(:,1),mu(:,2)); figure; 
%   pdfplot(fu)
%
% See also chitwo2lc_sorm, dirsp2chitwo


% References: U. Machado, I. Rychlik (2003) "Wave statistics in nonlinear
%                                                 sea", Extremes, 6, pp. 125--146.           
%             Butler, R., Machado, U. Rychlik, I. (2003) "Distribution of wave crests in non-
%             linear random sea - application of saddlepoint methods", presented at ISOPE 2003.
%             Hagberg, O. (2005) PhD - thesis, Dept. of Math. Statistics, Univ. of Lund. 
%   Calls: kchitwo
% Revised by I.R 14.04.05
% Revised by I.R 24.11.03
%
%------------------------------------------------------------------------------------



if nargin<1 
  error('problem is not defined');
end

n=length(gamma);
if nargin<2 
  beta=zeros(size(gamma));
end
variance=2*(sum(gamma.^2));
if nargin<3 
  S12=zeros(n);
end
if nargin<4 
  S22=eye(n);
end
if nargin<5 || ds0<=0
  ds0=0.1;
end

  ucrt=5*sqrt(sum(beta.^2)+2*variance);
 

% ====================================================================================


K=[];
DK=[];
DDK=[];
DK111=[];

DK1122=[];
DK122=[];

DDKt=[];

DK2222=[];
R0=[];
S=[];

s=0;
j=0;
flag=0;

ds1=ds0;

% In these formulas dt must be equal to ds

%==>positive side
K10old=0;
while flag~=1
 
  j=j+1;
  
  if (j>400)
      flag=1;
      disp('step dt is too small ')
%      error('step dt is too small - stop')
  end 
% ds=0.2*ds1;
 ds=0.05*ds1;
 dt=ds;
	 
 [K0]=kchitwo(s,0,gamma,beta,S12,S22);
 [K1]=kchitwo(s-ds,0,gamma,beta,S12,S22);
 [K2]=kchitwo(s+ds,0,gamma,beta,S12,S22);   
 [K3]=kchitwo(s-2*ds,0,gamma,beta,S12,S22);
 [K4]=kchitwo(s+2*ds,0,gamma,beta,S12,S22);
 
 [K1t]=kchitwo(s,-dt,gamma,beta,S12,S22);
 [K2t]=kchitwo(s,dt,gamma,beta,S12,S22);
 
 [K11]=kchitwo(s-ds,-dt,gamma,beta,S12,S22);
 [K22]=kchitwo(s+ds,dt,gamma,beta,S12,S22);
 
 [K12t]=kchitwo(s-ds,dt,gamma,beta,S12,S22);
 [K21t]=kchitwo(s+ds,-dt,gamma,beta,S12,S22);
 
 [K2tt]=kchitwo(s,2*dt,gamma,beta,S12,S22);
 [K1tt]=kchitwo(s,-2*dt,gamma,beta,S12,S22);

 K10=(K2-K1)/(2*ds);
 K20=(K1+K2-2*K0)/(ds*ds);
 
 K30=(K4-2*K2+2*K1-K3)/(2*ds*ds*ds);
 K40=(K4-4*K2+6*K0-4*K1+K3)/(ds*ds*ds*ds);

 K12=(K22-2*K2+K21t-K12t+2*K1-K11)/(2*ds*ds*ds);
 
 K02=(K1t+K2t-2*K0)/(dt*dt);
 K04=(K2tt-4*K2t+6*K0-4*K1t+K1tt)/(dt*dt*dt*dt);
 K22=(K22-2*K2t+K12t-2*K2+4*K0-2*K1+K21t-2*K1t+K11)/(ds*ds*ds*ds);  
%[j K10];
dds=abs(K10-K10old);
K10old=K10;
if abs(dds)>0.075
    ds1=ds1/2;
end
 if (j~=1) && ((s*K10-K0)>25 || 0>(s*K10-K0) || abs(K10)>ucrt || K20<0)

%if (j~=1) & (abs(K0-s*K10)>20)
%	if (j~=1) & (abs(K0-s*K10)>10) 
  flag=1;
else

 S=[S; s]; 
 s=s+ds1/2; 
 
 K=[K; K0];  
 R0=[R0; (K40/K20^2/8-5*K30*K30/K20^3/24)];
 DK=[DK; K10];
 DDK=[DDK; K20];
 
 DK111=[DK111; K30];
 DDKt=[DDKt; K02];
 DK1122=[DK1122; K22]; 
 DK122=[DK122; K12]; 
 DK2222=[DK2222; K04];
 end	
end	 

%=======================================================

s=-ds0/2;
j=0;
flag=0;
ds1=ds0;
K10old=DK(1);
%==>negative side
ds=0.25*ds1;
dt=ds;

while flag~=1
 j=j+1;

 
  if (j>400)
      flag=1;
      disp('step dt is too small ')
%      error('step dt is too small - stop')
  end 
 [K0]=kchitwo(s,0,gamma,beta,S12,S22);
 [K1]=kchitwo(s-ds,0,gamma,beta,S12,S22);
 [K2]=kchitwo(s+ds,0,gamma,beta,S12,S22);
 [K3]=kchitwo(s-2*ds,0,gamma,beta,S12,S22);
 [K4]=kchitwo(s+2*ds,0,gamma,beta,S12,S22);
  
 [K1t]=kchitwo(s,-dt,gamma,beta,S12,S22);
 [K2t]=kchitwo(s,dt,gamma,beta,S12,S22); 
 
 [K11]=kchitwo(s-ds,-dt,gamma,beta,S12,S22);
 [K22]=kchitwo(s+ds,dt,gamma,beta,S12,S22);

 [K12t]=kchitwo(s-ds,dt,gamma,beta,S12,S22);
 [K21t]=kchitwo(s+ds,-dt,gamma,beta,S12,S22);
 
 [K2tt]=kchitwo(s,2*dt,gamma,beta,S12,S22);
 [K1tt]=kchitwo(s,-2*dt,gamma,beta,S12,S22);
 
 K10=(K2-K1)/(2*ds);
 K20=(K1+K2-2*K0)/(ds*ds);
 
 K30=(K4-2*K2+2*K1-K3)/(2*ds*ds*ds); 
 K40=(K4-4*K2+6*K0-4*K1+K3)/(ds*ds*ds*ds);

 K12=(K22-2*K2+K21t-K12t+2*K1-K11)/(2*ds*ds*ds);
 
 K02=(K1t+K2t-2*K0)/(dt*dt);
 K04=(K2tt-4*K2t+6*K0-4*K1t+K1tt)/(dt*dt*dt*dt);
 K22=(K22-2*K2t+K12t-2*K2+4*K0-2*K1+K21t-2*K1t+K11)/(ds*ds*ds*ds); 
% [j K10];
 if (j~=1) && ((s*K10-K0)>25 || 0>(s*K10-K0) || abs(K10)>ucrt || K20<0)
%     ((s*K10-K0)>25 | 0>(s*K10-K0) )  & abs(K10)>ucrt  & K20<0
%		 if (j~=1) & (abs(K0-s*K10)>10) 
    flag=1;
else
dds=abs(K10-K10old);
K10old=K10;
if dds>0.1
    ds1=ds1/2;
end
 
 S=[s;S];
 s=s-ds1/2;
 ds=0.25*ds1;
 dt=ds;
 K=[K0;K];
 R0=[(K40/K20^2/8-5*K30*K30/K20^3/24); R0];

 DK=[K10;DK];
 DDK=[K20;DDK];
 
 DK111=[K30;DK111];
 DDKt=[K02; DDKt];
 DK1122=[K22; DK1122]; 
 DK122=[K12; DK122]; 
 DK2222=[K04; DK2222];
 
end
end	 

% ====================================================================================
% ====================================================================================
% marginal density

f=1/sqrt(2*pi)*exp(K-S.*DK)./sqrt(DDK);
wh=sign(S).*sqrt(2*(S.*DK-K));
%f1=pdfnorm(wh)./sqrt(DDK);

[area]=trapz(DK,f);
beta=sqrt(2*abs(S.*DK-K)-2*log(area));
%area=1;
%fu=[DK f/area]; 
fu=createpdf;
  Htxt = 'Density of X(0) non-central Chi^2 process';
  xtxt = 'X(0)';
fu.title=Htxt;
fu.labx{1}=xtxt;
fu.x{1}=DK;
fu.f=f/area;   

%fu1=[DK f1/area];
F1=cdfnorm(wh)-pdfnorm(wh).*((1./S)./sqrt(DDK)-1./wh);
Fu=[DK F1];

% ====================================================================================
% 1-TERM SADDLEPONT APRROXIMATION (GAUSSIAN APPR)
% ====================================================================================

muG=(f/area).*sqrt(DDKt)/sqrt(2*pi);
cc=sqrt(DDKt)/sqrt(2*pi);
%figure(3)
%plot(DK,muG,'k:o')

% ====================================================================================
% 2-TERMS SADDLEPONT APRROXIMATION (GAUSSIAN APPR)
% ====================================================================================

% equation 22
eq22=-1/4*(DK1122.*DDK-DK111.*DK122)./((DDK.^2).*DDKt);

%mu_2=muG.*(1+eq22);

%figure(3)
%hold on
%plot(DK,mu_2,'k:*')

% ====================================================================================
% 3-TERMS SADDLEPONT APRROXIMATION (GAUSSIAN APPR)
% ====================================================================================

% equation 23
eq23=(DK2222.*DDK-3*DK122.^2)./(DDK.*(DDKt.^2));

mu1=muG.*(1+eq22);
mu=muG.*(1+R0).*(1+eq22-(1/24)*eq23);
mu2=muG.*(1+eq22-(1/24)*eq23.*(1-eq22));
%figure(3)
%hold on
%plot(DK,mu,'k:d')

%mu1=DDKt;
%mu2=(DK2222.*DDK-3*DK122.^2)./DDK;

mu=[DK mu];
muG=[DK muG];
mu1=[DK mu1];
mu2=[DK mu2];





% ====================================================================================
