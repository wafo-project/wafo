function Sk=time2spa(S,k,k2,g,rate)
%TIME2SPA Transform of spectrum from frequency to wave no. (used in spec2spec)
%
% CALL:  Sk = time2spa(S,k,k2,g)
%
%   Sk = wave number spectrum, 1- or 2D depending on input
%   S  = spectrum struct
%  k,k2= wave number lag of output (default depending on input spectrum)
%   g  = constant of gravity (default see gravity)
% rate = defining the ratio between the wave numbers and angular
%        frequencies. (default 2)
%
% Transforms a frequency or directional spectrum to 
% a wave number spectrum.
% Input 1D gives output 1D, input 2D gives output 2D
% Rotation, transformation etc is preserved
%
% See also definitions, spec2spec

% Tested on Matlab 5.3
% History:  
%   revised pab feb 2007
%   - added physical rotation of spectrum
%   revised by så 03.10.2005: corrected bug when input k=dk
%   revised by es 29.11.1999: norm. in both x and y (length(g) 1 OR 2)
%   by es 99.08.17

error(nargchk(1,5,nargin))
if nargin<5|isempty(rate),
  rate = 2;
end
rate = max(round(abs(rate)),1);

if nargin<3
  if isfield(S,'g')
    g=S.g;
  else
    g=gravity;
  end
end
if strcmpi(S.type(end-2:end),'k1d')|strcmpi(S.type(end-2:end),'k2d')
  error('Spectrum already in space domain')
end
Sk = S;
if isfield(S,'f')
  w=2*pi*S.f(:);
  Sf=S.S/2/pi;
  Sk=rmfield(Sk,'f');
else
  w=S.w(:);
  Sf=S.S;
  Sk=rmfield(Sk,'w');
end
if isfield(Sk,'v') 
  %Sk=rmfield(Sk,{'v','phi'});
  % enc.velocity  and direction make no sense for wave no spectrum
  if Sk.v~=0
    error('enc.velocity  and direction make no sense for wave no spectrum')
  end
  Sk = rmfield(Sk,{'v'});
end
Sk.g=g;
nw = rate*length(w);
if strcmpi(S.type,'freq')|strcmpi(S.type,'enc')
  if nargin<2|isempty(k)
    k=linspace(0,w2k(w(end),0,S.h,g(1)),nw)';
  else
    if length(k)==1 % then interpret k as dk (step length)
      k=((0:length(w)-1)*k)';
    end
  end

  wk = k2w(k,0,Sk.h,g(1));
  in = (wk>min(w))&(wk<=w(end));
  Gk = zeros(size(w));
  Gk(in) = interp1(w,Sf,wk(in),'pchip');
  Sk.S   = zeros(size(k'));
  Sk.S   =(Gk.*dwdk(wk,k,Sk.h,g(1)))';
  Sk.type='k1D';
  Sk.k=k';
else % directional spectrum
  if S.h<inf
    % TODO % Make transformation from S(w,theta) -> S(k,k2) when h<inf
    error('WAFO:time2spa','This transformation for finite depth is not available yet')
  end
  nw = ceil(nw/2);
  if nargin<2|(isempty(k) & isempty(k2)) % no arg-in for wave-numbers
    k=linspace(0,w(end)^2/g(1),nw);
    if g(end)>0
      k2=linspace(-w(end)^2/g(end), w(end)^2/g(end), 2*nw-1)';
    else
      k2=linspace(-w(end)^2/g(1), w(end)^2/g(1), 2*nw-1)';
    end     
  else % Some arg-in for k
    if nargin<3|isempty(k2) % no arg-in for k2
      k2=k;
    end
    if length(k)==1 % k scalar, ie dk
      k=(0:(nw-1))*k;
    end
    if length(k2)==1 % ditto for k2
      k2 = (-(nw-1):(nw-1))'*k2;
    end
  end
  % k,k2 vectors. Make matrices
  k=k(:)'; % make k row vector
  if (S.theta(1)<-pi/2)||(S.theta(end)>pi/2);
    %thetas on the left half plane => k's<0
    k=[-fliplr(k(:,2:end)) k];
  end
  k2    = k2(:); % make k2 column vector
  %  K1    = ones(size(k2))*k;
  %  K2    = k2*ones(size(k));
  K1    = k(ones(length(k2),1),:);
  K2    = k2(:,ones(length(k),1));
  Sk.k  = k;
  Sk.k2 = k2;

  % Matrices K1 and K2 are equidistant
  W  = sqrt(sqrt(g(1)^2*K1.^2+g(end)^2*K2.^2)); % Non-equidistant W derived from (K1,K2)
  TH = atan2(g(end)*K2,g(1)*K1); % => -pi < TH <= pi, Non-eq.dist. TH ----"----
  Gk = zeros(size(W));
  in = (W>w(1))&(W<=w(end))&(TH>S.theta(1))&(TH<=S.theta(end));
  if 0,
    Gk(in) = interp2(w,S.theta,Sf,W(in),TH(in),'*linear');
  else
    % physically rotate grid!
    Gk(in) = interp2(w,S.theta,Sf,W(in),mod(TH(in)+S.phi+pi,2*pi)-pi,'*linear');
    Sk.phi = 0;
  end
  Sk.S   = zeros(size(Gk));
  Sk.S(W>0) = Gk(W>0)./abs(W(W>0)).^3/2*g(1)*g(end);  % Wave-number spectrum
  Sk        = rmfield(Sk,'theta');
  Sk.type  = 'k2D';
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dw=dwdk(w,k,h,g)

dw = zeros(size(w));
pw = w>0;
if h==inf % or k*h>15
  dw(pw)=g/2./w(pw);
else
  tanhkh = tanh(k(pw)*h);
  dw(pw)=(g*tanhkh+g*h*k(pw).*(1-tanhkh.^2))./w(pw)/2;
end
return









