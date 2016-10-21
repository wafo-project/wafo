function S1 = jonswap(w1,sdata,plotflag)
%JONSWAP Calculates (and plots) a JONSWAP spectral density
%
% CALL:  S = jonswap(w,sdata,plotflag); 
%        S = jonswap(wc,sdata,plotflag); 
%
%   S     = a struct containing the spectral density. See datastructures
%   w     = angular frequency        (default linspace(0,wc,257))
%   wc    = angular cutoff frequency (default 33/Tp)
%   sdata = [Hm0 Tp gamma sa sb A], where
%          Hm0   = significant wave height (default 7 (m))
%          Tp    = peak period (default 11 (sec))
%          gamma = peakedness factor determines the concentraton
%                  of the spectrum on the peak frequency,  1 <= gamma <= 7. 
%                  (default depending on Hm0, Tp, see below)
%          sa,sb = spectral width parameters (default 0.07 0.09)
%          Ag    = normalization factor used when gamma>1, (default -1)
%                  Ag<0  : Ag calculated by integration so that int S dw =Hm0^2/16
%                  Ag==0 : Ag = (1+1.00*log(gamma)^1.16)/gamma  
%                  Ag>0  : Ag = Ag
%       plotflag = 0, do not plot the spectrum (default).
%                  1, plot the spectrum.       
%
% JONSWAP return a spectrum object defined as
%
%         S(w) = A * Gf * G0 * wn^(-N)*exp(-N/(M*wn^M))
%    where 
%         G0  = Normalizing factor related to Bretschneider form
%         A   = Ag * (Hm0/4)^2 / wp     (Normalization factor)
%         Gf  = j^exp(-.5*((wn-1)/s)^2) (Peak enhancement factor) 
%         wn  = w/wp       
%         wp  = angular peak frequency
%         s   = sa      for wn <= 1 
%               sb      for 1  <  wn 
%         j   = gamma,     (j=1, => Bretschneider spectrum) 
%
%
%  The JONSWAP spectrum is assumed to be especially suitable for the North Sea, 
%  and does not represent a fully developed sea. It is a reasonable model for
%  wind generated sea when 3.6*sqrt(Hm0) < Tp < 5*sqrt(Hm0) 
%  A standard value for gamma is 3.3. However, a more correct approach is 
%  to relate gamma to Hm0:
%        D = 0.036-0.0056*Tp/sqrt(Hm0);
%        gamma = exp(3.484*(1-0.1975*D*Tp^4/(Hm0^2)));
%  This parameterization is based on qualitative considerations of deep water
%  wave data from the North Sea, see Torsethaugen et. al. (1984)
%  Here gamma is limited to 1..7.
%
%  The relation between the peak period and mean zero-upcrossing period 
%  may be approximated by
%         Tz = Tp/(1.30301-0.01698*gamma+0.12102/gamma)
%
%  For zero values, NaN's or parameters not specified in DATA the
%  default values are used. 
%
% Example:  % Bretschneider spectrum Hm0=7, Tp=11
%      S = jonswap([],[0 0 1])
%
% See also  pmspec, torsethaugen, simpson

% References:
% Torsethaugen et al. (1984)
% Characteristica for extreme Sea States on the Norwegian continental shelf. 
% Report No. STF60 A84123. Norwegian Hydrodyn. Lab., Trondheim
%
% Hasselmann et al. (1973)
% Measurements of Wind-Wave Growth and Swell Decay during the Joint
% North Sea Project (JONSWAP). 
% Ergansungsheft, Reihe A(8), Nr. 12, Deutschen Hydrografischen Zeitschrift.


% Tested on: matlab 6.0, 5.3
% History:
% Revised pab april 2007
% - replaced code with call to mkjonswap
% revised pab June 2005
% -fixed a bug in help header: the jonswap range is now correct
% revised pab 11jan2004
% - replaced code with call to getjonswappeakedness  
% revised jr 22.08.2001
% - correction in formula for S(w) in help: j^exp(-.5*((w-wp)/s*wp)^2) 
%   (the first minus sign added)
% revised pab 01.04.2001 
% - added wc to input
% revised jr 30.10 2000
%   - changed 'data' to 'sdata' in the function call
% revised pab 20.09.2000 
%   - changed default w: made it dependent on Tp
% revised es 25.05.00 
%   - revision of help text  
% revised pab 16.02.2000
%  - fixed a bug for sa,sb and the automatic calculation of gamma.
%  - added sa, sb, A=alpha to data input 
%  - restricted values of gamma to 1..7
% revised by pab 01.12.99
%  added gamma to data input
% revised by pab 11.08.99
% changed so that parameters are only dependent on the 
% seastate parameters Hm0 and Tp.
% also checks if Hm0 and Tp are reasonable.


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

%monitor=0; 

if nargin<3||isempty(plotflag),  plotflag=0;end

Hm0=7;Tp=11; gam=0; sa=0.07; sb=0.09; A=-1;% default values
data2=[Hm0 Tp gam sa sb A];
nd2=length(data2);
if (nargin>1) && ~isempty(sdata), 
  nd=length(sdata); 
  ind=find(~isnan(sdata(1:min(nd,nd2))));
  if any(ind) % replace default values with those from input data
    data2(ind)=sdata(ind);
  end
end
if (nd2>0) && (data2(1)>0),
  Hm0 = data2(1);
end
if (nd2>1) && (data2(2)>0),
  Tp = data2(2);
end
if (nd2>2) && (data2(3)>=1) && (data2(3)<=7), 
  gam = data2(3);
end
if (nd2>3) && (data2(4)>0),
  sa = data2(4);
end
if (nd2>4) && (data2(5)>0), 
  sb = data2(5);
end
if (nd2>5) ,
  A = data2(6);
end



w = [];
if nargin<1||isempty(w1), 
  wc = 33/Tp;
elseif length(w1)==1,
  wc = w1; 
else
  w = w1 ;
end
nw = 257;
if isempty(w), 
  w = linspace(0,wc,nw).';
end

options = mkjonswap('defaults');
options.Hm0 = Hm0;
options.Tp = Tp;
if gam>=1
  options.gamma = gam;
end
options.sigmaA = sa;
options.sigmaB = sb;
if A>0
  options.Ag = A;
elseif A==0
  options.method = 'parametric';
end
jspec = mkjonswap(options);
options = jspec('options');
gam = options.gamma;

S1   = createspec;
S1.S = jspec(w(:));
S1.w = w(:);
S1.norm =0; % The spectrum is not normalized
S1.note =['JONSWAP, Hm0 = ' num2str(Hm0)  ', Tp = ' num2str(Tp), ', gamma = ' num2str(gam)];




if plotflag
  plotspec(S1,plotflag)
end
