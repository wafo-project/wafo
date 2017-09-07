function [D,options] = mkspreading(type,varargin)
%MKSPREADING Function handle to Directional spreading function.
%
% CALL:  D = mkspreading(type,options);
%        D = mkspreading(type,Property1,value1,Property2,value2,....);
%
%      D    = function handle to a Directional spreading function returning
%            S = D(theta,w,wc) where S is a Nt X Nw matrix with the principal 
%           direction always along the x-axis. 
%      type = type of spreading function, see options below (default 'cos2s')
%         'cos2s'  : cos-2s spreading    N(S)*[cos((theta-theta0)/2)]^(2*S)  (0 < S) 
%         'box'    : Box-car spreading   N(A)*I( -A < theta-theta0 < A)      (0 < A < pi)
%         'mises'  : von Mises spreading N(K)*exp(K*cos(theta-theta0))       (0 < K)
%         'poisson': Poisson spreading   N(X)/(1-2*X*cos(theta-theta0)+X^2)  (0 < X < 1)
%         'sech2'  : sech-2 spreading    N(B)*sech(B*(theta-theta0))^2       (0 < B)
%         'wnormal': Wrapped Normal      
%                  [1 + 2*sum exp(-(n*D1)^2/2)*cos(n*(theta-theta0))]/(2*pi) (0 < D1)
%         (N(.) = normalization factor)       
%         (the first letter is enough for unique identification)  
%   options = Options structure with the fields:
%     .theta0   = function handle, inline object, matrix or a scalar defining 
%                 average direction given in radians at every angular frequency.
%                (length 1 or length == length(wn)) (default 0)
%     .method   = Defines function used for direcional spreading parameter:
%                 0, 'none'     : S(wn) = spa, frequency independent
%                 1, 'mitsuyasu': S(wn) frequency dependent  (default)
%                 2, 'donelan'  : B(wn) frequency dependent 
%                 3, 'banner'   : B(wn) frequency dependent 
%       S(wn) = spa *(wn)^ma,  : wnlo <= wn < wnc
%             = spb *(wn)^mb,  : wnc  <= wn < wnup
%             = 0              : wnup <= wn 
%       B(wn) = S(wn)          : wnlo <= wn < wnup
%             = spb*wnup^mb    : wnup <= wn,         method = 2
%             = sc*F(wn)       : wnup <= wn ,        method = 3
%               where F(wn) = 10^(-0.4+0.8393*exp(-0.567*log(wn^2))) and sc is
%               scalefactor to make the spreading funtion continous.
%      .wnlimits = [wnlo wnc wnup] limits used in the function defining the
%                 directional spreading parameter. 
%                 wnc is the normalized cutover frequency   (default [0 1 inf])
%      .sp       = [spa,spb] maximum spread parameters      (default [15 15])
%      .m        = [ma,mb] shape parameters                 (default [5 -2.5])
%
%  
% MKSPREADING Return function handle to a Directional spreading function.
% Here the S- or B-parameter, of the COS-2S and SECH-2 spreading function,
% respectively, is used as a measure of spread. All the parameters of the 
% other distributions are related to this parameter through the first Fourier
% coefficient, R1, of the directional distribution as follows: 
%         R1 = S/(S+1) or S = R1/(1-R1).
% where 
%         Box-car spreading  : R1 = sin(A)/A
%         Von Mises spreading: R1 = besseli(1,K)/besseli(0,K), 
%         Poisson spreading  : R1 = X
%         sech-2 spreading   : R1 = pi/(2*B*sinh(pi/(2*B))
%         Wrapped Normal     : R1 = exp(-D1^2/2)
%
% A value of S = 15 corresponds to 
%           'box'    : A=0.62,          'sech2'  : B=0.89      
%           'mises'  : K=8.3,           'poisson': X=0.94  
%           'wnormal': D=0.36
%
% The COS2S is the most frequently used spreading in engineering practice.
% Apart from the current meter/pressure cell data in WADIC all
% instruments seem to support the 'cos2s' distribution for heavier sea
% states, (Krogstad and Barstow, 1999). For medium sea states
% a spreading function between COS2S and POISSON seem appropriate,
% while POISSON seems appropriate for swell.
%   For the COS2S Mitsuyasu et al. parameterized SPa = SPb =
% 11.5*(U10/Cp) where Cp = g/wp is the deep water phase speed at wp and
% U10 the wind speed at reference height 10m. Hasselman et al. (1980)
% parameterized  mb = -2.33-1.45*(U10/Cp-1.17).
% Mitsuyasu et al. (1975) showed that SP for wind waves varies from 
% 5 to 30 being a function of dimensionless wind speed.
% However, Goda and Suzuki (1975) proposed SP = 10 for wind waves, SP = 25
% for swell with short decay distance and SP = 75 for long decay distance.
% Compared to experiments Krogstad et al. (1998) found that ma = 5 +/- eps and
% that -1< mb < -3.5. 
% Values given in the litterature:    [spa  spb   ma   mb      wlim(1:3)  ]
%   (Mitsuyasu: spa == spb)  (cos-2s) [15   15    5    -2.5  0    1    3  ]
%   (Hasselman: spa ~= spb)  (cos-2s) [6.97 9.77  4.06 -2.3  0    1.05 3  ]
%   (Banner   : spa ~= spb)  (sech2)  [2.61 2.28  1.3  -1.3  0.56 0.95 1.6] 
%
% Examples :
%  options = mkspreading('defaults'); % get default options.
%  options.sp(1) = 10;                % Set spa = 10
%  D = mkspreading('cos2s',options);
%  w = linspace(0,3,257);
%  theta = linspace(-pi,pi,129);
%  contour(D(theta,w))
%
%        % Make frequency dependent direction
%   options.theta0 = inline('pi/6*w');
%   D2 = mkspreading('cos2s',options);
%  contour(D2(theta,w))
%
%  % Plot all spreading functions
% alltypes = {'cos2s','box','mises','poisson','sech2','wnormal'};
% for ix =1:length( alltypes)
%   D3 = mkspreading(alltypes{ix},options);
%   figure(ix)
%   contour(D3(theta,w)),title(alltypes{ix})
% end
%
% close all
% 
% See also  mkdspec, plotspec, spec2spec
 
% References
%  Krogstad, H.E. and Barstow, S.F. (1999)
%  "Directional Distributions in Ocean Wave Spectra"
%  Proceedings of the 9th ISOPE Conference, Vol III, pp. 79-86
%
%  Goda, Y. (1999)
%  "Numerical simulation of ocean waves for statistical analysis"
%  Marine Tech. Soc. Journal, Vol. 33, No. 3, pp 5--14 
%
%  Banner, M.L. (1990)
%  "Equilibrium spectra of wind waves."
%  J. Phys. Ocean, Vol 20, pp 966--984
%
% Donelan M.A., Hamilton J, Hui W.H. (1985)
% "Directional spectra of wind generated waves."
% Phil. Trans. Royal Soc. London, Vol A315, pp 387--407
%
% Hasselmann D, Dunckel M, Ewing JA (1980)
% "Directional spectra observed during JONSWAP."
%  J. Phys. Ocean, Vol.10, pp 1264--1280
%
%  Mitsuyasu, H, et al. (1975)
%  "Observation of the directional spectrum of ocean waves using a
%  coverleaf buoy."
%  J. Physical Oceanography, Vol.5, No.4, pp 750--760

% Some of this might be included in help header:
% cos-2s:
% NB! The generally strong frequency dependence in directional spread
% makes it questionable to run load tests of ships and structures with a
% directional spread independent of frequency (Krogstad and Barstow, 1999).

% Parameterization of B
%    def = 2 Donelan et al freq. parametrization for 'sech2'
%    def = 3 Banner freq. parametrization for 'sech2'
%    (spa ~= spb)  (sech-2)  [2.61 2.28 1.3  -1.3  0.56 0.95 1.6] 


% Tested on: Matlab 7.0
% History:
% revised pab jan 2007
%  - renamed from spreading to mkspreading
%  - the function now return a function handle to the actual spreading function.
%  - removed wc, the cut over frequency -> input is now assumed as normalized frequency, w/wc.
% revised pab 17.06.2001
% - added wrapped normal spreading
% revised pab 6 April 2001
%  - added fourier2distpar
%  - Fixed the normalization of sech2 spreading
% revised by PAB and IR 1 April 2001: Introducing the azymuth as a
% standard parameter in order to avoid rotations of the directions
% theta. The x-axis is always pointing into the principal direction
% as defined in the spreading function D(omega,theta). The actual
% principal direction is defined by means of field D.phi.
% revised es 06.06.2000, commented away: if ((ma==0) & (mb==0)), ...,
%                    hoping that the check is unnecessary
% revised pab 13.06.2000
%  - fixed a serious bug: made sure -pi<= th-th0 <=pi
% revised pab 16.02.2000
%  -fixed default value for Hasselman parametrization 
% revised pab 02.02.2000
%   - Nt or th may be specified + check on th
%   - added frequency dependence for sech-2
%   - th0 as separate input
%   - updated header info
%   - changed check for nargins
%   - added the possibility of nan's in data vector
% Revised by jr 2000.01.25
% - changed check of nargins
% - frequency dependence only for cos-2s
% - updated information
% By es, jr 1999.11.25


%error(nargchk(0,inf,nargin));
narginchk(0,inf)
options = struct('theta0',0,'method','mitsuyasu',...
  'sp',[15 15],'m',[5 -2.5],'wnlimits',[0 1 inf]);

methodlist = {'none','mitsuyasu','donelan','banner'};

if (nargin==1) && strcmpi(type,'defaults')
  D = options;
  return
end
options = parseoptions(options,varargin{:});


if length(options.wnlimits)~=3
  error('WAFO:MKSPREADING','option.wnlimits must have 3 elements');
end

if length(options.m)~=2
  error('WAFO:MKSPREADING','option.m must have 3 elements');
end

% if isempty(options.wc) 
%   options.method = 'none';
% else
if isnumeric(options.method);
  options.method = methodlist{options.method+1};
end

methodNr = find(options.method(1)=='nmdb');
if isempty(methodNr)
  error('WAFO:MKSPREADING','Unknown method')
end

options.method = methodlist{methodNr};
switch options.method
  case {'donelan','banner'}
    if  ~isfinite(options.wnlimits(3))
      error('WAFO:MKSPREADING','option.wnlimits(3) must be finite for method>1')
    end  
end

if nargin<1||isempty(type), 
  type = 'cos2s'; 
end
switch type(1)
  case 'c'
    D = @(varargin)cos2s(options, varargin{:});
  case 'b'
    D = @(varargin)box(options,varargin{:});
  case 'm'
    D = @(varargin)mises(options,varargin{:});
  case 'p' 
    D = @(varargin)poisson(options,varargin{:});
  case 's'
    D = @(varargin)sech2(options,varargin{:});
  case 'w'
    D = @(varargin)wnormal(options,varargin{:});
  otherwise
    error('WAFO:MKSPREADING','Unknown type of spreading function');
end

end % mkspreading function

%% Nested subfunctions
function s =  getspread(wn, options) 
% GETSPREAD Return spread parameter, S, of COS2S function 
%
% CALL S =  getspread(wn, options) 
%
% S       = spread parameter of COS2S functions
% wn       = normalized frequencies.
% options = input options structure
  switch options.method,
  case {0,'none','n'} %isempty(wn)
    % no frequency dependent spreading,
    % but possible frequency dependent direction
    s = options.sp(1);
  otherwise
    %wn   = w./options.wc; % normalized frequency
    wlim = options.wnlimits;
    spa  = options.sp(1);
    spb  = options.sp(2);
    ma   = options.m(1);
    mb   = options.m(2);
    
    k = find(wlim(3)<wn);
    % Mitsuyasu et. al and Hasselman et. al parametrization   of
    % frequency dependent spreading
    numZeros = find(wn<=wlim(1),1,'last');
    s = [ zeros(numZeros,1) ; ...
      spa*(wn((wlim(1)<wn) & (wn<=wlim(2)))).^ma ;...
      spb*(wn((wlim(2)<wn) & (wn<=wlim(3)))).^mb ;...
      zeros(length(k),1) ].';
    if any(k)
      switch options.method
        case {2,'donelan','d'}
          % Donelan et. al. parametrization for B in SECH-2
          s(k) = spb*(wlim(3)).^mb;
          
         % Convert to S-paramater in COS-2S distribution
          r1 = r1ofsech2(s);
          s  = r1/(1-r1);
          
        case {3,'banner','b'}
          % Banner parametrization  for B in SECH-2
         s3m   = spb*(wlim(3)).^mb;
         s3p   = donelan(wlim(3));
         scale = s3m/s3p; % Scale so that parametrization will be continous
         s(k)  = scale*donelan(wn(k));
         
         % Convert to S-paramater in COS-2S distribution
          r1 = r1ofsech2(s);
          s  = r1/(1-r1);
         
      end
    end
  end
end % get spread


  function [s_par,TH,phi0,Nt] = inputchk(options, theta,w,wc, type)
    % INPUTCHK
    %
    %  CALL [s_par,TH,phi0,Nt] = inputchk(options, theta,w,wc, type)
    %

    % Default values
    %~~~~~~~~~~~~~~~
%   Nt = 101;
%   Nw = 257;
% 
%   if nargin<1 || length(theta)<2
%     if (nargin>0) && (length(theta)==1) && round(theta)==theta,
%       Nt = abs(theta);
%     end
%     theta = linspace(-pi,pi,Nt).';
%   elseif abs(theta(end)-theta(1))<2*pi-eps | abs(theta(end)-theta(1))>2*pi+eps
%     error('WAFO:mkspreading','theta must cover all angles -pi -> pi')
%   else
%     Nt = length(theta);
%     theta      = theta(:);
%   end
%   if Nt<40,
%     warning('WAFO:mkspreading','Number of angles is less than 40. Spreading too sparsely sampled!'),
%   end
   if nargin<4 || isempty(wc),
     wc = 1; % Peak frequency
   end
  theta = theta(:);
  Nt = length(theta);


  if nargin<3||isempty(w),
    w    = 1 ; %linspace(0,6,Nw).';
  else
    w   = w(:)/wc;
  end

  % Make sure theta is from -pi to pi
  if 1,
    phi0 = 0;
    theta = mod(theta+pi,2*pi)-pi;
  else
    phi0       = theta(1)+pi;
    theta      = theta-phi0;
    theta(end) = pi;
    phi0       = mod(phi0,2*pi);
  end  
 
 
  switch class(options.theta0)
    case {'function_handle','inline'}
      th0 = options.theta0(w).';
      if size(th0,1)>1
        th0 = th0.';
      end
    otherwise
      if isnumeric(options.theta0)
        th0 = options.theta0(:).';
      else
        error('WAFO:MKSPREADING','options.theta0 contains an unsupported datatype');
      end
  end
  
  Nt0 = length(th0);
  Nw  = length(w);
  isFreqDepDir = Nt0==Nw;
  if isFreqDepDir,
    % frequency dependent spreading and/or
    % frequency dependent direction
    % make sure -pi<=TH<pi
    TH = mod(theta(:,ones(1,Nw))-th0(ones(Nt,1),:)+pi,2*pi)-pi; 
  elseif Nt0~=1,
    error('WAFO:MKSPREADING','The length of theta0 must equal to 1 or the length of w');
  else

    TH = mod(theta-th0+pi,2*pi)-pi; % make sure -pi<=TH<pi
    if options.method(1)~='n', % frequency dependent spreading
      TH = TH(:,ones(1,Nw));
    end
  end
  
  % Get spreading parameter
  %~~~~~~~~~~~~~~~~
  s = getspread(w,options);
  if any(s<0), 
    error('WAFO:MKSPREADING','The COS2S spread parameter, S(w), value must be larger than 0');
  end
  
  if type(1)=='c'
    % cos2s 
    s_par = s;
  else
    r1 = abs(s./(s+1)); % First Fourier coefficient of the directional spreading function.
    
    % Find distribution parameter from first Fourier coefficient.
    s_par = fourier2distpar(r1,type);
  end
end % function inputchk


function [D, phi0] = cos2s(options, theta,w,wc)
%COS2S COS2S spreading function
% 
% CALL  [D, phi0] = cos2s(options, theta,w,wc);
%       [D, phi0] = cos2s(Nt,w); 
% 
%   D    = matrix of directonal spreading function, COS2S, defined as
%
%         cos2s(theta,w) =  N(S)*[cos((theta-theta0)/2)]^(2*S)  (0 < S) 
%
%         where N() is a normalization factor and S is the spreading parameter 
%         possibly dependent on w. The principal direction of D is always along 
%         the x-axis. size Nt X Nw. 
%   phi0 = Parameter defining the actual principal direction of D.
%  theta = vector of angles given in radians. Lenght Nt
%  w     = vector of angular frequencies given rad/s. Length Nw.
%
    if strncmpi(theta,'options',1)
      D = options;
      return
    end

    if nargin<3
      w = [];
    end
    if nargin<4
      wc = [];
    end

  
    [s,TH,phi0,Nt] = inputchk(options, theta, w, wc, 'cos2s');

    
    if options.method(1)=='n'
      S=s;
    else
      S = s(ones(Nt,1),:);
    end
    D = (exp(gammaln(S+1)-gammaln(S+1/2))/(2*sqrt(pi))).*cos(TH/2).^(2*S);  
  end % function cos2s



  function [D,phi0] = poisson(options, theta,w,wc)
    %POISSON POISSON spreading function
    %
    % CALL  [D, phi0] = poisson(options, theta,w,wc);
    %
    %   D    = matrix of directonal spreading function, POISSON, defined as
    %
    %         poisson(theta,w) =   N(X)/(1-2*X*cos(theta-theta0)+X^2)  (0 < X < 1)
    %
    %         where N() is a normalization factor and X is the spreading parameter
    %         possibly dependent on w. The principal direction of D is always along
    %         the x-axis. size Nt X Nw.
    %   phi0 = Parameter defining the actual principal direction of D.
    %
    %
    
    if strncmpi(theta,'options',1)
      D = options;
      return
    end
    if nargin<3
      w = [];
    end
    if nargin<4
      wc = [];
    end

    
    [r1,TH,phi0,Nt] = inputchk(options, theta,w,wc, 'poisson');
    
    if options.method(1)=='n'
      X = r1;
    else
      X = r1(ones(Nt,1),:);
    end
    if any(X>=1),
      error('WAFO:MKSPREADING:POISSON','POISSON spreading: X value must be less than 1');
    end
    D    = (1-X.^2)./(1-(2*cos(TH)-X).*X)/(2*pi);
  end % function poissong

  function [D,phi0] = wnormal(options, theta,w,wc)
    %WNORMAL Wrapped Normal spreading function
    %
    % CALL  [D, phi0] = wnormal(options, theta,w,wc);
    %
    %   D    = matrix of directonal spreading function, WNORMAL, defined as
    %
    %         wnormal(theta,w) = N(D1)*[1 + 2*sum exp(-(n*D1)^2/2)*cos(n*(theta-theta0))]  (0 < D1)
    %
    %         where N() is a normalization factor and D1 is the spreading parameter
    %         possibly dependent on w. The principal direction of D is always along
    %         the x-axis. size Nt X Nw.
    %   phi0 = Parameter defining the actual principal direction of D.
    %
    
    if strncmpi(theta,'options',1)
      D = options;
      return
    end
    if nargin<3
      w = [];
    end
    if nargin<4
      wc = [];
    end

    [par,TH,phi0,Nt] = inputchk(options, theta,w,wc, 'wnormal');

    D1 = par.^2/2;
    if options.method(1)~='n',
      D1 = D1(ones(Nt-1,1),:);     
    end
    
    ix  = (1:Nt-1).';
    ix2 = ix.^2;
    Nd2 = size(D1,2);
    Fcof = [ ones(1,Nd2)/2;...
      exp(-ix2(:,ones(1,Nd2)).*D1)     ]/pi;
    Pcor = [ ones(1,Nd2 ); exp(sqrt(-1)*ix*TH(1,:))]; % correction term to get
    % the correct integration limits
    Fcof = Fcof(1:Nt,:).*conj(Pcor);
    D = real(fft(Fcof));
    D(D<0)=0;
  end % function wrapped normal.

  function [D,phi0] = sech2(options, theta,w,wc)
    %SECH2 Sech^2 spreading function
    %
    % CALL  [D, phi0] = sech2(options, theta,w,wc);
    %
    %   D    = matrix of directonal spreading function, SECH2, defined as
    %
    %         sech2(theta,w) = N(B)*0.5*B*sech(B*(theta-theta0))^2 (0 < B)
    %
    %         where N() is a normalization factor and X is the spreading parameter
    %         possibly dependent on w. The principal direction of D is always along
    %         the x-axis. size Nt X Nw.
    %   phi0 = Parameter defining the actual principal direction of D.
    %
    
    if strncmpi(theta,'options',1)
      D = options;
      return
    end
    if nargin<3
      w = [];
    end
    if nargin<4
      wc = [];
    end

    [par,TH,phi0,Nt] = inputchk(options, theta,w,wc, 'sech2');
    NB = tanh(pi*par); % Normalization factor.
    NB(NB==0) = 1;     % Avoid division by zero
    if options.method(1)=='n'
      B=par;
    else
      B = par(ones(Nt,1),:);
      NB = NB(ones(Nt,1),:);
    end
    D = 0.5*B.*sech(B.*TH).^2./NB;
  end % function sech-2

 function [D,phi0] = mises(options, theta,w,wc)
   %MISES Mises spreading function
   %
   % CALL  [D, phi0] = mises(options, theta,w,wc);
   %
   %   D    = matrix of directonal spreading function, MISES, defined as
   %
   %         mises(theta,w) =  N(K)*exp(K*cos(theta-theta0))       (0 < K)
   %
   %         where N() is a normalization factor and K is the spreading parameter
   %         possibly dependent on w. The principal direction of D is always along
   %         the x-axis. size Nt X Nw.
   %   phi0 = Parameter defining the actual principal direction of D.
   %
  
   if strncmpi(theta,'options',1)
     D = options;
     return
   end
   if nargin<3
     w = [];
   end
   if nargin<4
     wc = [];
   end

   [par,TH,phi0,Nt] = inputchk(options, theta,w,wc, 'mises');
   if options.method(1)=='n'
     K = par;
   else
     K = par(ones(Nt,1),:);
   end
   D = exp(K.*(cos(TH)-1))./(2*pi*besseli(0,K,1));
   
   
 end % function mises

 function [D,phi0] = box(options, theta,w,wc)
   %BOX Box car spreading function
   %
   % CALL  [D, phi0] = box(options, theta,w,wc);
   %
   %   D    = matrix of directonal spreading function, SECH2, defined as
   %
   %         box(theta,w) =  N(A)*I( -A < theta-theta0 < A)      (0 < A < pi)
   %
   %         where N() is a normalization factor and A is the spreading parameter
   %         possibly dependent on w. The principal direction of D is always along
   %         the x-axis. size Nt X Nw.
   %   phi0 = Parameter defining the actual principal direction of D.
   %
   % BOX Return the 
   
if strncmpi(theta,'options',1)
  D = options;
  return
end
if nargin<3
  w = [];
end
if nargin<4
  wc = [];
end

  [par,TH,phi0,Nt] = inputchk(options, theta,w,wc, 'box');
  if any(par>pi), 
    error('WAFO:MKSPREADING:BOX','BOX-CAR spreading: The A value must be less than pi');
  end
  if options.method(1)=='n'
    A = par; 
  else
    A = par(ones(Nt,1),:);
  end
  D=(-A<=TH & TH<=A)./(2.*A); 
 end % function box
 %% End of Nested functions


%% Local sub functions

function x=fourier2distpar(r1,type)
%FOURIER2DISTPAR Fourier coefficients to distribution parameter  
%
% CALL x = fourier2distpar(r1,type)
%
% x    = distribution parameter
% r1   = corresponding fourier coefficient.
% type = string defining spreading function
%        'box'  
%        'mises'
%        'poisson'
%        'sech2'
%        'wnormal'
% 
% The S-parameter of the COS-2S spreading function is used as a measure of
% spread in MKSPREADING. All the parameters of the other distributions are
% related to this S-parameter through the first Fourier coefficient, R1, of the
% directional distribution as follows: 
%         R1 = S/(S+1) or S = R1/(1-R1).
% where 
%         Box-car spreading  : R1 = sin(A)/A
%         Von Mises spreading: R1 = besseli(1,K)/besseli(0,K), 
%         Poisson spreading  : R1 = X
%         sech-2 spreading   : R1 = pi/(2*B*sinh(pi/(2*B))
%         Wrapped Normal     : R1 = exp(-D1^2/2)


x = zeros(size(r1));

switch type(1),
  case 'p', %Poisson
     x = r1;
  case 'w', % wrapped normal
     x(r1<=0) = inf;
     ind = find(r1>0); % avoid log of 0
     x(ind) = sqrt(-2*log(r1(ind)));
    %D1 = sqrt(-2*log(r1));
    if any(x<0|~isreal(x)),
      error('WAFO:MKSPREADING:WNORMAL''WRAPPED NORMAL spreading: D value must be greater than 0.'),
    end
    
  case 's', % sech2
    B0 = [linspace(eps,5,513) linspace(5.0001,100)]';
    B0 = interp1(r1ofsech2(B0),B0,r1(:));
    k1 = find(isnan(B0));
    if any(k1),
      B0(k1) = 0;
      B0(k1) = max(B0);
    end  
    ix0 = find(r1~=0);

    for ix=ix0,
      x(ix) = abs(fzero(@(x) 0.5*pi./(sinh(.5*pi./x))-x.*r1(ix),B0(ix)));
    end
  case 'm', % mises
    x0 = [linspace(0,10,513) linspace(10.00001,100)]';
    x0 = interp1(besseli(1,x0,1)./besseli(0,x0,1),x0,r1(:));
    k1 = find(isnan(x0));
    if any(k1),
      x0(k1)= 0;
      x0(k1)=max(x0);
    end
  
    
    for ix=1:length(r1),
      x(ix) = fzero(@(x) besseli(1,x,1)./besseli(0,x,1)-r1(ix),x0(ix));
    end
    
  case 'b', % Box-char
    % Initial guess
    x0 = linspace(0,pi+0.1,1025).';
    x0 = interp1(sinc(x0/pi),x0,r1(:));
    x  = x0.';
    
    % Newton-Raphson 
    dx = ones(size(r1));
    iy=0;
    max_count = 100;
    ix = find(x);
    while (any(ix) && iy<max_count)
      xi=x(ix);
      dx(ix) = (sin(xi) -xi.*r1(ix)) ./(cos(xi)-r1(ix));
      xi = xi - dx(ix);
      % Make sure that the current guess is larger than zero and less than pi
      %x(ix) = xi + 0.1*(dx(ix) - 9*xi).*(xi<=0) + 0.38*(dx(ix)-6.2*xi +6.2).*(xi>pi) ;
      % Make sure that the current guess is larger than zero.
      x(ix) = xi + 0.5*(dx(ix) - xi).* (xi<=0);
      
      iy=iy+1;
      ix=find((abs(dx) > sqrt(eps)*abs(x))  &  abs(dx) > sqrt(eps)); 
      %disp(['Iteration ',num2str(iy),...
      %'  Number of points left:  ' num2str(length(ix)) ]),
    end
    if iy == max_count, 
      warning('WAFO:MKSPREADING:FOURIER2DISTPAR','Newton raphson method did not converge.');
      str = 'The last step was:  ';
      outstr = sprintf([str,'%13.8f'],dx(ix));
      fprintf(outstr);
    end
    x(x==0)=eps; % Avoid division by zero  
end
  %semilogy(x,abs(x(:)-x0))
  %plot(x,x(:)-x0)
end % fourier2distpar function


function s = donelan(wn)
%DONELAN High frequency decay of B of sech2 paramater
s =  10.^(-0.4+0.8393*exp(-0.567*log(wn.^2)));
end


function y = r1ofsech2(B)
%R1OFSECH2   Computes R1 = pi./(2*B.*sinh(pi./(2*B)));

x = 2.*B;

y           = zeros(size(x));
isXLarge    = (realmax<=x);
y(isXLarge) = 1;

k    = (realmin<x & x<100);
xk   = pi./x(k);
y(k) = xk./sinh(xk);

k2    = (100<=x & x<realmax);
if any(k2)
  xk2   = pi./x(k2);
  y(k2) = -2*xk2./(exp(xk2).*expm1(-2*xk2));
end
end % function r1ofsech2

function y = sinc(x)
%SINC Sin(pi*x)/(pi*x) function.
%
%  CALL:  y = sinc(x);
%
%    y = sin(pi*x)/(pi*x)    if x ~= 0
%      = 1                   if x == 0
%
% Example:
%  x =linspace(-5,5)';
%  plot(x,sinc(x))
%
% See also  sin

%Tested on: Matlab 5.3
%History:
% by pab 05.04.2001

y = ones(size(x));
k = find(x~=0);
if any(k);
  xk = pi*x(k);
  y(k) = sin(xk)./(xk);
end



end % sinc function

