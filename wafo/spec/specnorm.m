function [S,mn4,m0,m2,m4,m1]=specnorm(spectrum,m0,m2,plotflag)
%SPECNORM Normalize a spectral density such that m0=m2=1
%
% CALL:  [Sn,mn4,m0,m2,m4,m1]=specnorm(S,m0,m2,plotflag);
%
%        Sn  = the normalized spectrum struct
%        mn4 = the fourth spectral moment of Sn.
%        mi  = the i'th spectral moment of S.
%        S   = the original spectrum struct
%   plotflag = 0, do not plot the normalized spectrum (default).
%              1, plot the normalized spectrum.
%           
% Normalization performed such that
%    INT S(freq) dfreq = 1       INT freq^2  S(freq) dfreq = 1
% where integration limits are given by  freq  and  S(freq)  is the 
% spectral density; freq can be frequency or wave number.
% The normalization is defined by
% A=sqrt(m0/m2); B=1/A/m0; freq'=freq*A; S(freq')=S(freq)*B;
%
% If S is a directional spectrum then a normalized gravity (.g) is added
% to Sn, such that mxx normalizes to 1, as well as m0 and mtt.
% (See spec2mom for notation of moments)
%
% If S is complex-valued cross spectral density which has to be
% normalized, then m0, m2 (suitable spectral moments) should be given.
%
% Example: 
%   S = jonswap;
%   [Sn,mn4] = specnorm(S);
%   assert(spec2mom(Sn,2), [1,1], 1e-4)     % Should be equal to one!
% 
% See also  spec2mom, spec2spec 

% Tested on: Matlab 5.3
%
% History:
% revised pab
% -renamed from wnormspec -> specnorm
% revised ir 31.08.01, introducing normalizing  g for w-spectrum
% Revised by jr 20.04.2001
% - Condition on nargout related to calculation of mn4; minor change
% - Updated information, added example
% By es 28.09.1999

if  nargin==2 && ~isempty(m0)
  plotflag=m0;
end
if (nargin<4 && nargin ~=2)||isempty(plotflag)
  if nargout==0,
    plotflag=1;
  else
    plotflag=0;
  end
end

if strcmpi(spectrum.type(end-2:end),'k2d')
  intype=spectrum.type;
  spectrum=spec2spec(spectrum,'dir');
end
  
S=spectrum.S;   % size np x nf

if isfield(spectrum,'w')
  f=spectrum.w(:); % length nf
elseif isfield(spectrum,'k')
  f=spectrum.k(:);
else
  f=2*pi*spectrum.f(:);
  S=S/2/pi;
end

if (nargin<3)||isempty(m0)||isempty(m2)
  S1=abs(S);
  if strcmpi(spectrum.type(end-2:end),'dir')
    % integrate out theta:
    S2=trapz(spectrum.theta(:),...
	     S1.*(cos(spectrum.theta(:)*ones(size(f'))).^2),1);
   
    m20=trapz(f,f.^4/gravity^2.*S2.',1); % second order moment in x
    m02=trapz(f,f.^4/gravity^2.*trapz(spectrum.theta(:),S1.*...
	(sin(spectrum.theta(:)*ones(size(f'))).^2),1).',1); % second order moment in y
    S1=trapz(spectrum.theta(:),S1,1).'; % integrate out theta
    mom=trapz(f,[S1  f.*S1 f.^2.*S1 f.^4.*S1] ,1);
    m0=mom(1); m1=mom(2); m2=mom(3);  m4=mom(4);
  else
    S1=S1(:);
    mom=trapz(f,[S1  f.*S1 f.^2.*S1 f.^4.*S1] ,1);
    m0=mom(1); m1=mom(2); m2=mom(3);  m4=mom(4);
    if isfield(spectrum,'w')
      m02=m4/gravity^2;
      m20=m02;
    end
  end
end
 
SM0 = sqrt(m0); SM2=sqrt(m2);
A = SM0/SM2;
B = SM2/(SM0*m0);

f  = f*A;
S1 = S*B;

if (nargout >= 2) && nargin<3
  mn4 = m4*A^5*B;
elseif (nargout >= 2) && (nargin>3)
  mn4=trapz(f, f.^4.*S1);
end

S=spectrum;
S.S=S1;
if isfield(spectrum,'w')
  S.w=f;
elseif isfield(spectrum,'k')
  S.k=f;
else
  S.f=f/2/pi;
  S.S=S.S*2*pi;
end
if strcmpi(spectrum.type(end-2:end),'dir')
  S.g=[gravity*sqrt(m0*m20)/m2 gravity*sqrt(m0*m02)/m2];
  % normalization gravity in x and y, resp.
end
if isfield(spectrum,'w')
  S.g = gravity*sqrt(m0*m20)/m2;
end

S.norm=1;
S.date=datestr(now);
if exist('intype','var')
  S=spec2spec(S,intype);
end  
if plotflag
  plotspec(S);
end
