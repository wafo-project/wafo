function H = mktmaspec(varargin)
%MKTMASPEC Function handle to JONSWAP spectrum for finite water depth
%
% CALL:  H = mktmaspec(options,par1,val1,par2,val2,...) 
%
% H      = function handle to modified JONSWAP spectral density. 
%    options = structure with the fields:
%     .h     = water depth             (default 42 (m))
%     .Hm0   = significant wave height (default  7 (m))
%     .Tp    = peak period (default 11 (sec))
%     .gamma = peakedness factor determines the concentraton
%              of the spectrum on the peak frequency.
%              Usually in the range  1 <= gamma <= 7. 
%              default depending on Hm0, Tp, see getjonswappeakedness)
%     .sigmaA = spectral width parameter for w<wp (default 0.07)
%     .sigmaB = spectral width parameter for w<wp (default 0.09)
%     .Ag     = normalization factor used when gamma>1: 
%     .N      = scalar defining decay of high frequency part.   (default 5)
%     .M      = scalar defining spectral width around the peak. (default 4)
%     .method = String defining method used to estimate Ag when gamma>1
%              'integrate' : Ag = 1/gaussq(Gf.*ggamspec(wn,N,M),0,wnc) (default)
%              'parametric': Ag = (1+f1(N,M)*log(gamma)^f2(N,M))/gamma
%              'custom'    : Ag = Ag
%     .wnc    = wc/wp normalized cut off frequency used when calculating Ag 
%               by integration (default 6)
%
%
%  The evaluated spectrum is
%         S(w) = Sj(w)*phi(w,h)
%    where 
%         Sj  = jonswap spectrum
%         phi = modification due to water depth
%
%  The concept is based on a similarity law, and its validity is verified
%  through analysis of 3 data sets from: TEXEL, MARSEN projects (North
%  Sea) and ARSLOE project (Duck, North Carolina, USA). The data include
%  observations at water depths ranging from 6 m to 42 m.
%
% Example:  
%    options = mktmaspec('defaults'); % default options
%    options.h = 20;options.gamma = 1;
%    S = mktmaspec(options)   % Bretschneider spectrum Hm0=7, Tp=11
%    for h = [10 21 42]
%      S = mktmaspec('h',h);
%      fplot(S,[0 3]), hold on
%    end
% 
% See also  mkjonswap, phi1, mkbretschneider, mktorsethaugen

% References:
% Buows, E., Gunther, H., Rosenthal, W., and Vincent, C.L. (1985)
% 'Similarity of the wind wave spectrum in finite depth water: 1 spectral form.' 
%  J. Geophys. Res., Vol 90, No. C1, pp 975-986
%
% Hasselman et al. (1973)
% Measurements of Wind-Wave Growth and Swell Decay during the Joint
% North Sea Project (JONSWAP). 
% Ergansungsheft, Reihe A(8), Nr. 12, deutschen Hydrografischen
% Zeitschrift.

% Tested on: matlab 6.0, 5.2
% History:
% Revised pab Feb 2007
% -Function now return a function handle to spectrum
% -renamed from tmaspec -> mktmaspec
% revised jr 03.04.2001
% - added wc to input
% - updated information
% revised pab 22.09.2000, corrected water depth, i.e. changed 
%    S1.h=inf to S1.h=h
% revised es 25.05.00 - some modifications to help text    
% by pab 16.02.2000

 
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
%      especially important when studying velocities and accelerations.

error(nargchk(0,inf,nargin))


options = struct('h',42,'Hm0',7,'Tp',11,'gamma',[],...
  'sigmaA',0.07,'sigmaB',0.09,'Ag',[],'N',5,'M',4,...
  'method','integration','wnc',6,'chkseastate','on');

if (nargin==1) && ischar(varargin{1}) && strcmpi(varargin{1},'defaults')
  H = options;
  return
end
options = parseoptions(options,varargin{:});
Hj = mkjonswap(options); % Handle to jonswap function
H  = @tmaspec;

  function S = tmaspec(w)
  %TMASPEC Return actual TMA spec values  
  % 
  % CALL      S = tmaspec(w)
  %     options = tmaspec('options')
  %
 
    if strncmpi(w,'options',1)
      S = options;
    else
      S = Hj(w).*phi1(w,options.h);
    end
  end
end % function mktmaspec

function tmaspec()
%TMASPEC 
% This is a trick to get matlab help function work for nested functions!
%  help mktmaspec>tmaspec will work. 
end

