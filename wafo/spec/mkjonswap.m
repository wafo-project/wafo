function H = mkjonswap(varargin)
%MKJONSWAP Function handle to JONSWAP spectral density
%
% CALL:  H = mkjonswap(options,par1,val1,par2,val2,...)
%
%   H      = function handle to JONSWAP spectral density. 
%  options = structure with the fields:
%   .Hm0   = significant wave height (default 7 (m))
%   .Tp    = peak period (default 11 (sec))
%   .gamma = peakedness factor determines the concentraton
%            of the spectrum on the peak frequency.
%            Usually in the range  1 <= gamma <= 7. 
%            default depending on Hm0, Tp, see getjonswappeakedness)
%   .sigmaA = spectral width parameter for w<wp (default 0.07)
%   .sigmaB = spectral width parameter for w<wp (default 0.09)
%   .Ag     = normalization factor used when gamma>1: 
%   .N      = scalar defining decay of high frequency part.   (default 5)
%   .M      = scalar defining spectral width around the peak. (default 4)
%   .method = String defining method used to estimate Ag when gamma>1
%            'integrate' : Ag = 1/gaussq(Gf.*ggamspec(wn,N,M),0,wnc) (default)
%            'parametric': Ag = (1+f1(N,M)*log(gamma)^f2(N,M))/gamma
%            'custom'    : Ag = Ag
%   .wnc    = wc/wp normalized cut off frequency used when calculating Ag 
%             by integration (default 6)
%
% MKJONSWAP return a function handle to a JONSWAP spectrum defined as
%
%         S(w) = A * Gf * G0 * wn^(-N)*exp(-N/(M*wn^M))
%    where 
%         G0  = Normalizing factor related to Bretschneider form
%         A   = Ag * (Hm0/4)^2 / wp     (Normalization factor)
%         Gf  = j^exp(-.5*((wn-1)/s)^2) (Peak enhancement factor) 
%         wn  = w/wp       
%         wp  = angular peak frequency
%         s   = sigmaA      for wn <= 1 
%               sigmaB      for 1  <  wn 
%         j   = gamma,     (j=1, => Bretschneider spectrum) 
%
%  The JONSWAP spectrum is assumed to be especially suitable for the North Sea, 
%  and does not represent a fully developed sea. It is a reasonable model for
%  wind generated sea when the seastate is in the so called JONSWAP range, i.e., 
%    3.6*sqrt(Hm0) < Tp < 5*sqrt(Hm0) 
%
%  The relation between the peak period and mean zero-upcrossing period 
%  may be approximated by
%         Tz = Tp/(1.30301-0.01698*gamma+0.12102/gamma)
%
% Example:  % Bretschneider spectrum Hm0=7, Tp=11
%     options = mkjonswap('defaults');
%     assert(fieldnames(options), {'Hm0', 'Tp', 'gamma', 'sigmaA', 'sigmaB',...
%             'Ag', 'N', 'M', 'method', 'wnc', 'chkseastate'}')
%     assert(struct2cell(options), ...
%            {7,11,[],0.07,0.09,[],5,4,'integration', 6, 'on'}')
%     options.gamma = 1;
%     S = mkjonswap(options);
%  % or alternatively
%     S = mkjonswap('gamma', 1);
%  % Plot the spectrum by
%     fplot(S,[0,5]) 
%  % or alternatively
%     w = linspace(0,5);
%     plot(w, S(w))
% 
%   options2 = S('options'); % get options used
%
%   S2 = mkbretschneider(options);
%   [x,y] = fplot(S,[0,4]);
%   y2 = S2(x);
%   assert(y,y2,eps) %JONSWAP with gamma=1 equals Bretscneider!
%
%   close all
%
% See also  mkjonswap>jonswap, mktmaspec, mkbretschneider, mktorsethaugen, getjonswappeakedness

% References:
% Torsethaugen et al. (1984)
% Characteristica for extreme Sea States on the Norwegian continental shelf. 
% Report No. STF60 A84123. Norwegian Hydrodyn. Lab., Trondheim
%
% Hasselmann et al. (1973)
% Measurements of Wind-Wave Growth and Swell Decay during the Joint
% North Sea Project (JONSWAP). 
% Ergansungsheft, Reihe A(8), Nr. 12, Deutschen Hydrografischen Zeitschrift.


% Tested on: matlab 7, 8, 9
% History:
% By pab Jan 2007
% New enhanced implementation based on old jonswap and torsethaugen functions.



%NOTE: In order to calculate the short term statistics of the response,
%      it is extremely important that the resolution of the transfer
%      function is sufficiently good. In addition, the transfer function
%      must cover a sufficietn range of wave periods, especially in the
%      range where the wave spectrum contains most of its
%      energy. VIOLATION OF THIS MAY LEAD TO MEANINGLESS RESULTS FROM THE 
%      CALCULATIONS OF SHORT TERM STATISTICS. The highest wave period
%      should therefore be at least 2.5 to 3 times the highest peak
%      period in the transfer function. The lowest period should be selected 
%      so that the transfer function value is low. This low range is 
%      especially important when studying velocities a


%error(nargchk(0,inf,nargin))
narginchk(0,inf)

options = struct('Hm0',7,'Tp',11,'gamma',[],...
  'sigmaA',0.07,'sigmaB',0.09,'Ag',[],'N',5,'M',4,...
  'method','integration','wnc',6,'chkseastate','on');

if (nargin==1) && ischar(varargin{1}) && strcmpi(varargin{1},'defaults')
  H = options;
  return
end
options = parseoptions(options,varargin{:});
if isempty(options.gamma) || isnan(options.gamma) || options.gamma<1
  options.gamma = getjonswappeakedness(options.Hm0,options.Tp);
end
options = preCalculateAg(options);
if strcmpi(options.chkseastate,'on')
  chkseastate(options)
end


H = @(w)jonswap(options, w);
return
end % mkjonswap

%% Subfunctions

function S = jonswap(options, w)
%JONSWAP spectral density
% 
% CALL      S = jonswap(w)
%     options = jonswap('options')
%


  if strncmpi(w,'options',1)
    S = options;
  else
    if (options.Hm0>0)
    
      N   = options.N;
      M   = options.M;
      wp  = 2*pi/options.Tp;
      wn  = w/wp;
      Ag  = options.Ag;
      Hm0 = options.Hm0;
      Gf  = peakEnhancementFactor(wn,options);
      S   = ((Hm0/4)^2/wp.*Ag).*Gf.*ggamspec(wn,N,M);
    else
      S = zeros(size(w));
    end
  end
end % jonswap


function Gf = peakEnhancementFactor(wn,options)
%PEAKENHANCEMENTFACTOR 
%
% 

gam = options.gamma;
sb  = options.sigmaB;
sa  = options.sigmaA;


k    = (wn>1);
sab      = zeros(size(wn));
sab(k)   = sb;
sab(~k)  = sa;
wn(wn<0) = 0;

wnm12 = 0.5*((wn-1)./sab).^2;
Gf    = gam.^(exp(-wnm12));

end % peak enhancement factor

function options = preCalculateAg(options) 
% PRECALCULATEAG Precalculate normalization.

if (options.gamma==1)
   options.Ag = 1;
   options.method = 'parametric';
elseif ~isempty(options.Ag)
  options.method = 'custom';
  if options.Ag<=0
    error('WAFO:MKJONSWAP','Ag must be larger than 0!')
  end
elseif options.method(1)=='i'
  % normalizing by integration
   options.method = 'integration';
   if options.wnc<1
     error('WAFO:MKJONSWAP','Normalized cutoff frequency, wnc, must be larger than one!')
   end
   area = gaussq(@(x)localspec(x, options),0,1) + ...
          gaussq(@(x)localspec(x, options),1,max(options.wnc,0));
   % area = gaussq(@localspec,1,options.wc);
   options.Ag = 1/area;
elseif options.method(1)=='p'
  options.method = 'parametric';
  % Original normalization
   % NOTE: that  Hm0^2/16 generally is not equal to intS(w)dw
   %       with this definition of Ag if sa or sb are changed from the
   %       default values
   N = options.N;
   M = options.M;
   gammai = options.gamma;
   parametersOK = (3<=N && N<=50) &&  (2<=M && M <=9.5) && (1<= gammai && gammai<=20);
   if parametersOK
     f1NM = 4.1*(N-2*M^0.28+5.3)^(-1.45*M^0.1+0.96);
     f2NM = (2.2*M^(-3.3) + 0.57)*N^(-0.58*M^0.37+0.53)-1.04*M^(-1.9)+0.94;
     options.Ag = (1+ f1NM*log(gammai)^f2NM)/gammai;
     
%    elseif N == 5 && M == 4,
%      options.Ag = (1+1.0*log(gammai).^1.16)/gammai;
%      %options.Ag = (1-0.287*log(gammai));
%      options.normalizeMethod = 'Three';
%    elseif  N == 4 && M == 4,
%      options.Ag = (1+1.1*log(gammai).^1.19)/gammai;
   else
      error('WAFO:MKJONSWAP','Not knowing the normalization because N, M or peakedness parameter is out of bounds!')
   end
   if options.sigmaA~=0.07 || options.sigmaB~=0.09
     warning('WAFO:MKJONSWAP','Use integration to calculate Ag when sigmaA~=0.07 or sigmaB~=0.09')
   end
end
end % precalculateAg

function S = localspec(wn, options)
    Gf = peakEnhancementFactor(wn,options);
    S = Gf.*ggamspec(wn,options.N, options.M);
end % local spec

function chkseastate(options)
%CHKSEASTATE Check that seastate is OK
Tp  = options.Tp;
Hm0 = options.Hm0;
gam = options.gamma;

if Hm0<0
  error('WAFO:MKJONSWAP','Hm0 can not be negative!')
end

if Tp<=0
  error('WAFO:MKJONSWAP','Tp must be positve!')
end


if Hm0==0
  warning('WAFO:MKJONSWAP','Hm0 is zero!')
end

outsideJonswapRange  = Tp>5*sqrt(Hm0) || Tp<3.6*sqrt(Hm0);
if outsideJonswapRange
  warning('WAFO:MKJONSWAP','Hm0,Tp is outside the JONSWAP range \n The validity of the spectral density is questionable')
end
if gam<1 || 7 < gam
  warning('WAFO:MKJONSWAP','The peakedness factor, gamma, is possibly too large. \n The validity of the spectral density is questionable')
end
end % chkseastate