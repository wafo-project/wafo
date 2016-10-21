function [mu,fu,Fu,muG,mu1,mu2]= specq2lc(Spec,ds0,nb)
%SPECQ2LC Saddlepoint approximation of crossing intensity for quadratic sea.
%
%  CALL:  [mu,fu,Fu]= specq2lc(S,ds);
%  
%     mu   = three column vector containing: levels in the first column, 
%            the saddlepoint approximation of crossing intensity for quadratic sea  
%            in the second column and for Gaussian sea in the third column.
%            Note that the levels $u$ are not equally spaced. If the spacing 
%            in the grid is not fine enough, choose smaller value of parameter  ds. 
%      fu  = the pdf structure: f.x contains levels u, f.f is a two column matrix containing
%            the Danniel's-saddlepoint approximation of the density (pdf) of the quadratic 
%            sea in the first collon and pdf of the Gaussian sea in the second.  
%      Fu  = the two column matrix containing levels u and the Lugannani Rice's-saddlepoint 
%            approximation of the cdf of the quadratic sea.  
%     S    = spectral density structure, computed using WAFO. S.h contains the water depth.
%            (default: Jonswap with depth 5000 [m]).
%     ds   = parameter defining the levels: default value is ds=0.1.
%
% Example: 
%  S=jonswap; S.h=40;
%  [mu, fu, Fu]=specq2lc(S,0.05); semilogy(mu(:,1),mu(:,2:3))
%  pdfplot(fu)
%

% References: U. Machado, I. Rychlik (2003) "Wave statistics in non-linear random sea" 
%             Extremes, 6, pp. 125-146.
%             Butler, R., Machado, U. Rychlik, I. (2002): "Distribution of wave crests in non-
%             linear random sea - application of saddlepoint methods" by , presented at ISOPE 2003. 
% Calls: kwaves.m  which evaluates the cumulant generating function K(s,t) as defined in Eq. 13
% By U.M. 04.11.02
% Revised by I.R 22.11.02
%
%------------------------------------------------------------------------------------



if nargin<1 
  Spec=jonswap;
end

if nargin<2 || ds0<=0
  ds0=0.1;
end

if nargin<3
  nb=0;
end

h=Spec.h;
   if (h>5000)
       h=5000;
   end

%g=gravity;
%omega=Spec.w;
spec=sqrt(Spec.S);
spec=spec*spec';
dw=Spec.w(2,1)-Spec.w(1,1);
w=Spec.w;

if nb==1 % if narrow-band
 [i,j]=max(Spec.S);
 w=Spec.w(j)*ones(size(Spec.w));
end 

%=== Computation of the transfer functions Q, R and S

W=diag(-Spec.w);
Q=(eww(w,-w,h)+eww(w,w,h)).*spec*dw;
R=(eww(w,-w,h)-eww(w,w,h)).*spec*dw;
      
%===
Q=(Q+Q')/2;
R=(R+R')/2;
%S=Q*W-W*R;
%sigma=sqrt(Spec.S*dw);
%N=length(W);

%===
[P1c,Delta]=eig(Q);
P1=P1c';
% eigenvectors per line
%check: mesh(Q-P1'*Delta*P1)

[P2c,Gama]=eig(R);
P2=P2c';
%check: mesh(R-P2'*Delta*P2)

%==============================================
%              Delta P1 Gama P2
%
% ORDER from the small to the largest values
%==============================================
%
% small letters are vectors:

lambda=diag(real(Delta));
[Vl, Il]=sort(abs(lambda));
lambdaord(:,1)=lambda(Il(:,1),1);

%Vl: absolute values of eigenvalues in ascending order
%Il: positions where they were

%Deltaord=diag(lambdaord);
P1ord =P1(Il,:);
% for i=1:length(P1)
%  P1ord(i,:)=P1(Il(i,1),:);
% end
%====

gamda=diag(real(Gama));
[Vl, Ill]=sort(abs(gamda));
gamdaord(:,1)=gamda(Ill(:,1),1);

%Gamaord=diag(gamdaord);

P2ord = P2(Ill,:);

% for i=1:length(P2)
%  P2ord(i,:)=P2(Ill(i,1),:);
% end

%=== Computation of P now ordered
%P=[P1ord zeros(size(P1)) ;zeros(size(P2)) P2ord];

%==============================================
% REDUCE i.e. take away lines (replace mm lines by zeros)
%==============================================

variance=2*(sum(lambdaord.^2) + sum(gamdaord.^2));
varapprox=2*(cumsum(lambdaord.^2) + cumsum(gamdaord.^2))/variance;
mm=max(1,sum(varapprox<0.00001));

lambdatilde=lambdaord;
lambdatilde(1:mm,1)=zeros(mm,1);

gamdatilde=gamdaord;
gamdatilde(1:mm,1)=zeros(mm,1);

Deltatilde=diag(lambdatilde);
Gamatilde=diag(gamdatilde);

% position of the first element that is not zero
%np=mm+1;


% the same for both
nn=mm;
%npp=np;


% ====================================================================================

% for the computation r

sigma=sqrt(Spec.S*dw);
N=length(W);

% Here P1 is ordered and without zeros

rr1=P1ord*sigma;
rr2=P2ord*W*sigma;

r=[rr1;rr2];

% x1=zeros(mm,1);
% x2=zeros(N-mm,1);
% x3=zeros(nn,1);
% x4=zeros(N-nn,1);

x1=r(1:mm,1);
x2=r(mm+1:N,1);
x3=r(N+1:N+nn,1);
x4=r(N+nn+1:2*N,1);

sum1=sum(x1.^2)/2;
sum2=sum(x3.^2)/2;

Sn=Deltatilde*P1ord*W*P2ord'-P1ord*W*P2ord'*Gamatilde;
% ====================================================================================
%check: mesh(P1ord'*Sn*P2ord-S)

A22=Deltatilde(mm+1:N,mm+1:N);
A44=Gamatilde(nn+1:N,nn+1:N);
A24=Sn(mm+1:N,mm+1:N);
%A14=Sn(1:mm,mm+1:end);

C=zeros(2*mm,2*N-2*mm);
C(1:mm,N-mm+1:end)=Sn(1:mm,mm+1:end);
C(mm+1:end,1:N-mm)=Sn(mm+1:end,1:mm)';

CC=C'*C;


%%%%%%%%%%%%%%%%%%%%%%% CHECK of  accuracy of K %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nb==99
%  s1=-3;
%  s2=2;
%  t=[s1*sigma ; s2*W*sigma];
%  A=[s1*Q s2*S;  s2*S' s1*R];
%  IA=eye(size(A));
%  InvA=inv(IA-A);
%  IQ=inv(eye(size(Q))-s1*Q);
  %Note that we have  minus -S' giving the plus in the following formula.
%  r=s2*W*sigma+s2*S'*IQ*sigma*s1;
  % Here I am checking Lemma 6.
 % K0=0.5*s1^2*sigma'*IQ*sigma +0.5*r'*inv(eye(size(R))-s1*R-s2*S'*IQ*S*s2)*r-0.5*log(det(eye(size(Q))-s1*Q))-0.5*log(det(eye(size(R))-s1*R-s2*S'*IQ*S*s2));
 %   lamb = spec2mom(Spec);
%    [i,j]=max(Spec.S);
 %   wp=Spec.w(j);
 %   Const=wp*wp/(2*g);
%[ K0 0.5*t'*InvA*t-0.5*log(det(IA-A)) kwaves(s1,s2,A22,A44,A24,C',CC,sum1,sum2,x1,x2,x3,x4)   knb(s1,s2,lamb(1),lamb(2),lamb(3),Const)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% ====================================================================================


K=[];
DK=[];
DDK=[];
DK111=[];

DK1122=[];
DK122=[];

DDKt=[];

DK2222=[];

S=[];

s=0;
j=0;
flag=0;

ds=0.01;
dt=ds;

% In these formulas dt must be equal to ds

%==>positive side
while flag~=1
 
  j=j+1; 
  
	 
 [K0]=kwaves(s,0,A22,A44,A24,C',CC,sum1,sum2,x1,x2,x3,x4);
 [K1]=kwaves(s-ds,0,A22,A44,A24,C',CC,sum1,sum2,x1,x2,x3,x4);
 [K2]=kwaves(s+ds,0,A22,A44,A24,C',CC,sum1,sum2,x1,x2,x3,x4);   
 [K3]=kwaves(s-2*ds,0,A22,A44,A24,C',CC,sum1,sum2,x1,x2,x3,x4);
 [K4]=kwaves(s+2*ds,0,A22,A44,A24,C',CC,sum1,sum2,x1,x2,x3,x4);
 
 [K1t]=kwaves(s,-dt,A22,A44,A24,C',CC,sum1,sum2,x1,x2,x3,x4);
 [K2t]=kwaves(s,dt,A22,A44,A24,C',CC,sum1,sum2,x1,x2,x3,x4);
 
 [K11]=kwaves(s-ds,-dt,A22,A44,A24,C',CC,sum1,sum2,x1,x2,x3,x4);
 [K22]=kwaves(s+ds,dt,A22,A44,A24,C',CC,sum1,sum2,x1,x2,x3,x4);
 
 [K12t]=kwaves(s-ds,dt,A22,A44,A24,C',CC,sum1,sum2,x1,x2,x3,x4);
 [K21t]=kwaves(s+ds,-dt,A22,A44,A24,C',CC,sum1,sum2,x1,x2,x3,x4);
 
 [K2tt]=kwaves(s,2*dt,A22,A44,A24,C',CC,sum1,sum2,x1,x2,x3,x4);
 [K1tt]=kwaves(s,-2*dt,A22,A44,A24,C',CC,sum1,sum2,x1,x2,x3,x4);

 dK=(K2-K1)/(2*ds);
 ddK=(K1+K2-2*K0)/(ds*ds);
 
 dddK=(K4-2*K2+2*K1-K3)/(2*ds*ds*ds);
 d2sd2t=(K22-2*K2t+K12t-2*K2+4*K0-2*K1+K21t-2*K1t+K11)/(ds*ds*ds*ds); 
 dsd2t=(K22-2*K2+K21t-K12t+2*K1-K11)/(2*ds*ds*ds);
 
 ddKt=(K1t+K2t-2*K0)/(dt*dt);
 d4t=(K2tt-4*K2t+6*K0-4*K1t+K1tt)/(dt*dt*dt*dt);
 
 if (j~=1) && ((s*dK-K0)>25 || 0>(s*dK-K0) )
%if (j~=1) & (abs(K0-s*dK)>20)
%	if (j~=1) & (abs(K0-s*dK)>10) 
  flag=1;
else

 S=[S; s]; 
 s=s+ds0/2; 
 
 K=[K; K0];  
 DK=[DK; dK];
 DDK=[DDK; ddK];
 
 DK111=[DK111; dddK];
 DDKt=[DDKt; ddKt];
 DK1122=[DK1122; d2sd2t]; 
 DK122=[DK122; dsd2t]; 
 DK2222=[DK2222; d4t];
 end	
end	 

%=======================================================

s=-ds0/2;
j=0;
flag=0;

%==>negative side

while flag~=1
 j=j+1;
 [K0]=kwaves(s,0,A22,A44,A24,C',CC,sum1,sum2,x1,x2,x3,x4);
 [K1]=kwaves(s-ds,0,A22,A44,A24,C',CC,sum1,sum2,x1,x2,x3,x4);
 [K2]=kwaves(s+ds,0,A22,A44,A24,C',CC,sum1,sum2,x1,x2,x3,x4);
 [K3]=kwaves(s-2*ds,0,A22,A44,A24,C',CC,sum1,sum2,x1,x2,x3,x4);
 [K4]=kwaves(s+2*ds,0,A22,A44,A24,C',CC,sum1,sum2,x1,x2,x3,x4);
  
 [K1t]=kwaves(s,-dt,A22,A44,A24,C',CC,sum1,sum2,x1,x2,x3,x4);
 [K2t]=kwaves(s,dt,A22,A44,A24,C',CC,sum1,sum2,x1,x2,x3,x4); 
 
 [K11]=kwaves(s-ds,-dt,A22,A44,A24,C',CC,sum1,sum2,x1,x2,x3,x4);
 [K22]=kwaves(s+ds,dt,A22,A44,A24,C',CC,sum1,sum2,x1,x2,x3,x4);

 [K12t]=kwaves(s-ds,dt,A22,A44,A24,C',CC,sum1,sum2,x1,x2,x3,x4);
 [K21t]=kwaves(s+ds,-dt,A22,A44,A24,C',CC,sum1,sum2,x1,x2,x3,x4);
 
 [K2tt]=kwaves(s,2*dt,A22,A44,A24,C',CC,sum1,sum2,x1,x2,x3,x4);
 [K1tt]=kwaves(s,-2*dt,A22,A44,A24,C',CC,sum1,sum2,x1,x2,x3,x4);
 
 dK=(K2-K1)/(2*ds);
 ddK=(K1+K2-2*K0)/(ds*ds);
 
 dddK=(K4-2*K2+2*K1-K3)/(2*ds*ds*ds); 
 d2sd2t=(K22-2*K2t+K12t-2*K2+4*K0-2*K1+K21t-2*K1t+K11)/(ds*ds*ds*ds);
 dsd2t=(K22-2*K2+K21t-K12t+2*K1-K11)/(2*ds*ds*ds);
 
 ddKt=(K1t+K2t-2*K0)/(dt*dt);
 d4t=(K2tt-4*K2t+6*K0-4*K1t+K1tt)/(dt*dt*dt*dt);
 
 if (j~=1) && ((s*dK-K0)>25 || 0>(s*dK-K0) )
%		 if (j~=1) & (abs(K0-s*dK)>10)
    flag=1;
else
 
 S=[s;S];
 s=s-ds0/2;
 K=[K0;K];
 DK=[dK;DK];
 DDK=[ddK;DDK];
 
 DK111=[dddK;DK111];
 DDKt=[ddKt; DDKt];
 DK1122=[d2sd2t; DK1122]; 
 DK122=[dsd2t; DK122]; 
 DK2222=[d4t; DK2222];
 end
end	 

% ====================================================================================
% ====================================================================================
% marginal density

f=1/sqrt(2*pi)*exp(K-S.*DK)./sqrt(DDK);
%(S.*DK-K)
wh=sign(S).*sqrt(2*(S.*DK-K));
%f1=pdfnorm(wh)./sqrt(DDK);
L0=spec2mom(Spec,2);
f1=1/sqrt(2*pi)*exp(-DK.^2/2/L0(1))./sqrt(L0(1));

[area]=trapz(DK,f);
%fu=[DK f/area]; 
fu=createpdf;
  Htxt = 'Density of X(0) for guadratic sea';
  xtxt = 'X(0) [m]';
fu.title=Htxt;
fu.labx{1}=xtxt;
fu.x{1}=DK;
fu.f=[f/area f1];   

%fu1=[DK f1/area];
F1=cdfnorm(wh)-pdfnorm(wh).*((1./S)./sqrt(DDK)-1./wh);
Fu=[DK F1];

% ====================================================================================
% 1-TERM SADDLEPONT APRROXIMATION (GAUSSIAN APPR)
% ====================================================================================

muG=f.*sqrt(DDKt)/sqrt(2*pi)/area;
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
mu=muG.*(1+eq22-(1/24)*eq23);
mu2=muG.*(1+eq22-(1/24)*eq23.*(1-eq22));
%figure(3)
%hold on
%plot(DK,mu,'k:d')

%mu1=DDKt;
%mu2=(DK2222.*DDK-3*DK122.^2)./DDK;

mu=[DK mu 1/2/pi*exp(-DK.^2/2/L0(1))*sqrt(L0(2)/L0(1))];
muG=[DK muG];
mu1=[DK mu1];
mu2=[DK mu2];



return % main

% ====================================================================================

function [out]=eww(omega,omegat,h)
% EWW: Computes values of the quadratic transfer function E, given in Eq. 6. in
%             R. Butler, U. Machado, I. Rychlik (2002). 
%
%  CALL:  [Eww]= eww(w,wt,h);
%  
%     Eww  = a matrix with the quadratic transfer function E(w,wt). 
%     w,wt = two equally long vectors with angular frequencies.
%     h    = water depth (default 5000 [m]).
%
%   Example: S=Jonswap; w=[-flipud(S.w) ;S.w]; Eww=eww(w,w); mesh(w,w,Eww)

%  Bugs:   E(w,-w)=0. (This should be checked?)
%
%----------------------------------------------------------------------
% References: R. Butler, U. Machado, I. Rychlik (2002) Distribution of waves crests in nonlinear 
%             random sea - applications of saddlepoint methods. ISOPE 2003. 
% Calls: kwaves.m
% By U.M. 10.11.02
% Revised by I.R 22.11.02

g=gravity;
eps0=0.000001;

if (length(omega)~=length(omegat))
   error('error in input to eww')    
end

if nargin<3 || h<=0
  h=5000;
end

[w wt]=meshgrid(omega,omegat);
wpl=w+wt;


ind=find(abs(w.*wt.*wpl)<eps0);
wpl(ind)=1;
w(ind)=1;
wt(ind)=1;
	 
kw=w2k(w,[],h,g);
kwt=w2k(wt,[],h,g);

Etop=g*(kw.*kwt)./(w.*wt)-(w.^2+wt.^2+w.*wt)/(2*g)+...
	    	(g/2)*(w.*kwt.^2+wt.*kw.^2)./(w.*wt.*wpl);

Ebottom=1-g*(kw+kwt)./wpl.^2.*tanh((kw+kwt)*h);

out=Etop./Ebottom-(g*kw.*kwt)./(2*w.*wt)+(w.^2+wt.^2+w.*wt)/(2*g);

out(ind)=0;
return %eww

%-----------------------------------------------------------------------

function [K]=kwaves(s1,s2,A22,A44,A24,CT,CC,sum1,sum2,x1,x2,x3,x4)
% KWAVES: help function which computes the cumulant generating function 
%         formula (36) in Machado and Rychlik (2002) 
%         "Wave statistics in nonlinear random sea" to appear in Extremes.
%_____________________________________________________________________________ 
%   Note that s1,s2 must be constants (vectors are not allowed):


% tested on matlab 6.1
% By U Machado, last update 20.10.02
% revised IR 19.11.02, changed name and added some checks and comments.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We tried to keep notation of the Appendix in Machado&Rychlik
% CT means C',  CC=C'*C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% r=r(s_1,s_2) and B=B(s_1,s_2)   are defined in Eq. 37.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r=[s1*x2; s2*x4]+s2*CT*[s1*x1; s2*x3];

%r=[s1*x2; s2*x4]-s2*CT*[s1*x1; s2*x3];

%===

Aux=s2^2*CC;

I=eye(size(Aux));

B=I-([s1*A22 s2*A24; s2*A24' s1*A44])-Aux;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K=K(s_1,s_2) This is formulas (36) in Machado&Rychlik          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K=-1/2*log(det(B))+s1^2*sum1+s2^2*sum2+1/2*r'*inv(B)*r;

%[B rn]
%[K s1^2*sum1+s2^2*sum2 1/2*rn'*inv(B)*rn]
return % kwaves









