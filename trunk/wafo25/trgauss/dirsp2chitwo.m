function [gam,bet,S12,S22]= dirsp2chitwo(s,w,L0,L2,th,h,eps,dthdw)
% DIRSP2CHITWO Parameters in non-central CHI-TWO process for directional Stokes waves. 
%
%  CALL:  [gamma,beta,S12,S22]= dirsp2chitwo(s,w,L0,L2,th,h,eps,dwth);
%  
%     s,w,th = spectral density s(w,th), where w is a vector with angular frequencies
%              th wave directions (frequences and angles are equally spaced) Note that
%              usually the spectrum s is truncated to exclude low and high frequencies.
%              The degree of truncation is measured by using the spectral
%              moments L0, L2 for untruncated spectrum.
%      L0,L2 = Two first spectral moments of the linear spectrum. These can
%              be smaller than the corresponding spectral moments of  s. 
%          h = water depth (default: 5000 [m]).
%        eps = all eigenvalues which have absolutvalue below  eps*variance of the
%              CHI2 process are replaced by zero. 
%       dwth = spacing in w and th vectors dtwh=dw*dth
%
% gamma,beta,S12,S22 parameters in CHI-TWO model.
%
% See also chitwo2lc_sorm, chitwo2lc_sp


% References: U. Machado, I. Rychlik (2002) "Wave statistics in nonlinear sea" to
%             appear in Extremes.
%             Butler, R., Machado, U. Rychlik, I. (2002) "Distribution of wave crests in non-
%             linear random sea - application of saddlepoint methods" by , presented at ISOPE 2003. 
% Calls: ewwdir
%             By I.R 24.10.04
%
%------------------------------------------------------------------------------------



if nargin<2
  Spec=jonswap;
  s=Spec.S;
  w=Spec.w;
  th=zeros(size(w));
  h=Spec.h;
  %nb=0;
  dthdw=w(3)-w(2);
end
if nargin<4
    dthdw=w(3)-w(2);
    L0=dthdw*sum(s);
    L2=dthdw*sum(w.^2.*s);
end
if nargin<5
    th=zeros(size(w));
end
if nargin<6
    h=5000;
end
   if (h>5000)
       h=5000;
   end
if nargin<7
    eps=0.00001;
end

if nargin<8
    dth=th(3)-th(2);
    dthdw=(w(3)-w(2))*dth;
end
if (dthdw<0.0000000001)
    dthdw=(w(3)-w(2));
end

%g=gravity;
%omega=w;
spec=sqrt(s);
spec=spec*spec';

% if narrow-band
% [i,j]=max(Spec.S);
% w=Spec.w(j)*ones(size(Spec.w));
 

%=== Computation of the transfer functions Q, R and S

W=diag(-w);
Q=(ewwdir(w,th,-w,th,h)+ewwdir(w,th,w,th,h)).*spec*dthdw;
R=(ewwdir(w,th,-w,th,h)-ewwdir(w,th,w,th,h)).*spec*dthdw;


%===
Q=(Q+Q');
R=(R+R');
%S=Q*W-W*R;
%sigma=sqrt(s*dthdw);
N=length(W);


%===
[P1c,Delta]=eig(Q);
P1=P1c';
[P2c,Gama]=eig(R);
P2=P2c';



gam=[diag(Delta,0)' diag(Gama,0)']/2;
variance=2*(sum(gam.^2));


% computations of beta

sigma=sqrt(s*dthdw);
bet=[(P1*sigma)' zeros(1,N)];


Z=zeros(N);
SS12=[Z -W;W Z];
SS22=[W.^2 Z; Z W.^2];
PP=[P1 Z;Z P2];
S12=PP*SS12*PP';
S22=PP*SS22*PP';

%SS=eye(2*N);
%[bet*SS*bet' bet*S22*bet'];

%
% In this part of program we are removing the quadratic processes with
% negligable gamma_i coefficients.%

n=2*N;
if (eps>0)
 [gammasort indexgamma]=sort(abs(gam));

 gammasort=gam(indexgamma);
 betasort=bet(indexgamma);
 S12sort=S12(indexgamma,indexgamma);
 S22sort=S22(indexgamma,indexgamma);
% betasort*S22sort*betasort';
 varapprox=2*(cumsum(gammasort.^2))/variance;
 mm=sum(varapprox<eps);
 if (mm>0)
    beta1=betasort(1:mm);
    bet=[sqrt(sum(beta1.^2)) betasort(mm+1:end)];
    gam=[0 gammasort(mm+1:end)];
    S12=zeros(n-mm+1,n-mm+1);
    S22=zeros(n-mm+1,n-mm+1);
    S22(1,1)=beta1*S22sort(1:mm,1:mm)*beta1'/bet(1)^2;
    S22(1,2:end)=beta1*S22sort(1:mm,mm+1:end)/bet(1);
    S22(2:end,2:end)=S22sort(mm+1:end,mm+1:end);
    S12(2:end,2:end)=S12sort(mm+1:end,mm+1:end);
    S12(1,2:end)=beta1*S12sort(1:mm,mm+1:end)/bet(1);
    S12(2:end,1)=-S12(1,2:end)';
    S22(2:end,1)=S22(1,2:end)';
 end
end
dL0=L0-sum(s)*dthdw;
dL2=L2-(s'*w.^2)*dthdw;
bet(1)=sqrt(bet(1).^2+dL0);
if(dL2>0.001)
S22(1,1)=S22(1,1)+dL2/dL0;
end