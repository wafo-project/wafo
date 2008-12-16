function [Sn,mn,mom]=scalespec(So,m0n,m2n,plotflag)
%SCALESPEC  Scale spectral density so that the moments equals m0,m2.  
%
%  CALL: [Sn,mn,mo]=scalespec(So,m0n,m2n,plotflag);
%
%        Sn  = the scaled spectrum struct
%        mn  = [m0n m1n m2n m4n] the spectral moments of Sn.
%        mo  = [m0o m1o m2o m4o] the spectral moments of So.
%        So  = the original spectrum struct
%   plotflag = 0, do not plot the normalized spectrum (default).
%              1, plot the normalized spectrum.
%           
%  The Scaling is performed so that
%
%    INT Sn(freq) dfreq = m0n  INT freq^2  Sn(freq) dfreq = m2n
%
%  where  integration limits is given by freq and Sn(freq) is the 
%  spectral density. freq can be frequency or wave number.
%  Default values for m0n=m2n=1. The normalization is defined by:
%
%      freq'=freq*A; S(freq')=S(freq)*B;
%  where
%      A=sqrt(m0o/m2o)/sqrt(m0n/m2n); B=m0n/(A*m0o); 
%
%  If S is a directional spectrum then a normalized gravity (.g) is added
%  to Sn, such that mxx normalizes to 1, as well as m0 and mtt.
%  (See spec2mom for notation of moments)
%
% Example: Transform spectra from a model scale
%
%   Hm0 = 0.133; Tp = 1.36;
%   Sj = jonswap(linspace(0,125,1025),[Hm0,Tp,3]); 
%   ch = spec2char(Sj,{'Hm0','Tm02','Ss'});
%       % to the corresponding spectrum with Hm0=12 and Ss=ch(3)
%   Ss=ch(3);Tm02=ch(2);Hm0b=12;
%   m0n = (Hm0b/4)^2; 
%   m2n = 4*pi^2*m0n*Hm0/(Tm02^2*Hm0b); 
%   Sn = scalespec(Sj,m0n,m2n,1);
%   ch2 = spec2char(Sn,{'Hm0','Tm02','Ss'})
%
% See also  specnorm, spec2mom
 
% Tested on: Matlab 5.3
% History:
% revised pab jan2004  
% by pab 20.09.2000
  
% TODO % Scaling is not correct for directional spectra.
% TODO % Needs testing for directional spectra.


if nargin<3||isempty(m2n),m2n=1;end
if nargin<2||isempty(m0n),m0n=1;end   
if (nargin<4)||isempty(plotflag)
  if nargout==0,
    plotflag=1;
  else
    plotflag=0;
  end
end

if strcmpi(So.type(end-2:end),'k2d')
  intype=So.type;
  So=spec2spec(So,'dir');
end
  
S = So.S;   % size np x nf
ftype = freqtype(So);
switch ftype 
case 'w',  f=So.w(:);  
case 'k',  f=So.k(:);
otherwise
  f=2*pi*So.f(:);
  S=S/2/pi;
end
S1=abs(S);
if strcmp(So.type(end-2:end),'dir')
   S2=trapz(So.theta(:),S1.*(cos(So.theta(:)*ones(size(f'))).^2),1).'; % integrate out theta
   S1=trapz(So.theta(:),S1,1).'; % integrate out theta
   m20=trapz(f,f.^4/gravity^2.*S2,1);
else
   S1=S1(:);
end
mom = trapz(f,[S1  f.*S1 f.^2.*S1 f.^4.*S1] ,1);
m0 = mom(1); 
%m1 = mom(2); 
m2 = mom(3);  
%m4 = mom(4);

SM0 = sqrt(m0); 
SM2 = sqrt(m2); 
SM0n=sqrt(m0n); 
SM2n=sqrt(m2n); 
A = SM0*SM2n/(SM2*SM0n);
B = SM2*SM0n*m0n/(SM0*m0*SM2n);

f  = f*A;
S1 = S*B;

if (nargout > 1)
   mn=trapz(f,[S1  f.*S1 f.^2.*S1 f.^4.*S1] );
end

Sn=So;
Sn.S=S1;

switch ftype 
case 'w',  Sn.w=f; 
case 'k',  Sn.k=f; 
otherwise
   Sn.f=f/2/pi;
   Sn.S=S.S*2*pi;
end

if strcmp(So.type,'dir')
  Sn.g=gravity*sqrt(m0*m20)/m2;
end

Sn.norm=(m0n==1 & m2n==1);
if Sn.norm,
   Sn.note=[Sn.note ' Scaled'];
else	
   Hm0  = 4*sqrt(m0n);
   Tm02 = 2*pi*SM0n/SM2n;
   Sn.note=[Sn.note,' Scaled to Hm0 =' num2str(Hm0),' Tm02 =' num2str(Tm02)]; 
end
   
%Sn.date=strvcat(datestr(now),Sn.date);
Sn.date=datestr(now);
   
if exist('intype','var')
  Sn=spec2spec(Sn,intype);
end  
if plotflag
  plotspec(Sn)
end
 