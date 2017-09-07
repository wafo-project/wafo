function [H, Hw, Hs]=mktorsethaugen(varargin)
%MKTORSETHAUGEN Function handle to double peaked (swell + wind) spectrum 
%
%  CALL: [H, Hw, Hs] = mktorsethaugen(options); 
%
%   H, Hs, Hw = function handles to spectrum for total, swell and wind, respectively.  
%   options  =  structure with the fields:
%   .Hm0   = significant wave height (default 7 (m))
%   .Tp    = peak period (default 11 (sec))
%   .wnc    = wc/wp normalized cut off frequency used when calculating Ag 
%             by integration (default 6)
%   .method = String defining method used to estimate normalization factors, Ag,
%             in the the modified JONSWAP spectra when gamma>1
%            'integrate' : Ag = 1/gaussq(Gf.*ggamspec(wn,N,M),0,wnc)
%            'parametric': Ag = (1+f1(N,M)*log(gamma)^f2(N,M))/gamma
%
% MKTORSETHAUGEN return a function handle to a double peaked (swell + wind)
% spectrum modelled as  S(w) = Ss(w) + Sw(w) where Ss and Sw are modified
% JONSWAP spectrums for swell and wind peak, respectively.
% The energy is divided between the two peaks according
% to empirical parameters, which peak that is primary depends on parameters.   
% The empirical parameters are found for classes of Hm0 and Tp,
% originating from a dataset consisting of 20 000 spectra divided
% into 146 different classes of Hm0 and Tp. (Data measured at the
% Statfjord field in the North Sea in a period from 1980 to 1989.)
% The range of the measured  Hm0 and Tp for the dataset
% are from 0.5 to 11 meters and from 3.5 to 19 sec, respectively.
%
% Preliminary comparisons with spectra from other areas indicate that 
% some of the empirical parameters are dependent on geographical location.
% Thus the model must be used with care for other areas than the
% North Sea and sea states outside the area where measured data 
% are available.
%
% Example: 
%  S = mktorsethaugen('Hm0',6, 'Tp',8);
%  fplot(S,[0,5]);
%  options = mktorsethaugen('defaults');
%  assert(fieldnames(options), {'Hm0', 'Tp', 'method', 'wnc', 'chkseastate'}')
%  assert(struct2cell(options), {7,11,'integration', 6, 'on'}')
%  options = S('options');
%  assert(fieldnames(options), {'Hm0','Tp','method','wnc','chkseastate',...
%                               'Hwoptions','Hsoptions'}')
%  assert(fieldnames(options.Hwoptions), {'Hm0', 'Tp', 'gamma', 'sigmaA', 'sigmaB',...
%             'Ag', 'N', 'M', 'method', 'wnc', 'chkseastate'}')
%  assert(fieldnames(options.Hsoptions), {'Hm0', 'Tp', 'gamma', 'sigmaA', 'sigmaB',...
%             'Ag', 'N', 'M', 'method', 'wnc', 'chkseastate'}')
%
%  close all
%
% See also  mkjonswap, mkbretschneider.

% References 
%  Torsethaugen, K. (2004)
%  "Simplified Double Peak Spectral Model for Ocean Waves"
%  In Proc. 14th ISOPE
%
%  Torsethaugen, K. (1996)
%  Model for a doubly peaked wave spectrum 
%  Report No. STF22 A96204. SINTEF Civil and Environm. Engineering, Trondheim
%
%  Torsethaugen, K. (1994)
%  'Model for a doubly peaked spectrum. Lifetime and fatigue strength
%  estimation implications.' 
%  International Workshop on Floating Structures in Coastal zone,
%  Hiroshima, November 1994.
%
%  Torsethaugen, K. (1993)
%  'A two peak wave spectral model.'
%  In proceedings OMAE, Glasgow


% Tested on Matlab 5.3, 7. 8. 9
% History: 
% By pab Jan 2007
% New enhanced implementation based on old torsethaugen function from 1999.


%error(nargchk(0,inf,nargin));
narginchk(0,inf)
options = struct('Hm0',7,'Tp',11,'method','integration','wnc',6,'chkseastate','on');

if (nargin==1) && ischar(varargin{1}) && strcmpi(varargin{1},'defaults')
  H = options;
  return
end
options = parseoptions(options,varargin{:});
if strcmpi(options.chkseastate,'on')
  chkseastate(options)
end

[Hw, Hs] = mktspec(options);
H = @(w)torsethaugen(options, Hw, Hs, w);

return
end % mktorsethaugen

%% Subfunctions

function S = torsethaugen(options, Hw, Hs, w)
%TORSETHAUGEN spectral density
%
% CALL        S = torsethaugen(w)
%       options = torsethaugen('options')
%
  if strncmpi(w,'options',1)
    S = options;
    S.Hwoptions = Hw(w);
    S.Hsoptions = Hs(w);
  else
    S = Hw(w)+Hs(w);
  end
end


function chkseastate(options)
%
Hm0 = options.Hm0;
Tp = options.Tp;

if Hm0<0
  error('WAFO:MKTORSETHAUGEN','Hm0 can not be negative!');
end

if Tp<=0
  error('WAFO:MKTORSETHAUGEN','Tp must be positve!');
end


if Hm0==0
  warning('WAFO:MKTORSETHAUGEN','Hm0 is zero!');
end


if Hm0>11 || Hm0>max((Tp/3.6).^2, (Tp-2)*12/11)
  warning('WAFO:TORSETHAUGEN','Hm0 is outside the valid range.\n The validity of the spectral density is questionable');
end
if Tp>20||Tp<3
  warning('WAFO:TORSETHAUGEN','Tp is outside the valid range.\n The validity of the spectral density is questionable');
end



end % chk seastate
function [Hwind,Hswell] = mktspec(options)
  % MKTSPEC Return function handles to swell and wind part of Torsethaugen spectrum
  %
  % CALL:  [Hw,Hs] = mktspec(options)
  %
  %  Hw,Hs   = function handles to the wind- and swell-generated parts of
  %           the Torsethaugen spectrum, respectively
  monitor = 0;
  Hm0 = options.Hm0;
  Tp  = options.Tp;
  gravity1 = gravity; % m/s^2
  
  % The parameter values below are found comparing the
  % model to average measured spectra for the Statfjord Field
  % in the Northern North Sea.
  Af  = 6.6;   %m^(-1/3)*sec
  AL  = 2;     %sec/sqrt(m)
  Au  = 25;    %sec
  KG  = 35;
  KG0 = 3.5;
  KG1 = 1;     % m
  r   = 0.857; % 6/7
  K0  = 0.5;   %1/sqrt(m)
  K00 = 3.2;

  M0  = 4;
  B1  = 2;    %sec
  B2  = 0.7;
  B3  = 3.0;  %m
  S0  = 0.08; %m^2*s
  S1  = 3;    %m

  % Preliminary comparisons with spectra from other areas indicate that
  % the parameters on the line below can be dependent on geographical location
  A10 = 0.7; A1 = 0.5; A20 = 0.6; A2 = 0.3; A3 = 6;

  Tf = Af*(Hm0)^(1/3);
  Tl = AL*sqrt(Hm0);   % lower limit
  Tu = Au;             % upper limit
  
  %Non-dimensional scales
  % New call pab April 2005
  El = min(max((Tf-Tp)/(Tf-Tl),0),1); %wind sea
  Eu = min(max((Tp-Tf)/(Tu-Tf),0),1); %Swell



  if Tp<Tf, % Wind dominated seas
    % Primary peak (wind dominated)
    Nw  = K0*sqrt(Hm0)+K00;             % high frequency exponent
    Mw  = M0;                           % spectral width exponent
    Rpw = min((1-A10)*exp(-(El/A1)^2)+A10,1);
    Hpw = Rpw*Hm0;                      % significant waveheight wind
    Tpw = Tp;                           % primary peak period
    % peak enhancement factor
    gammaw = KG*(1+KG0*exp(-Hm0/KG1))*(2*pi/gravity1*Rpw*Hm0/(Tp^2))^r;
    gammaw = max(gammaw,1);
    % Secondary peak (swell)
    Ns  = Nw;                % high frequency exponent
    Ms  = Mw;                % spectral width exponent
    Rps = sqrt(1-Rpw^2);
    Hps = Rps*Hm0;           % significant waveheight swell
    Tps = Tf+B1;
    gammas = 1;
  
    if monitor
      if Rps > 0.1
        disp('     Spectrum for Wind dominated sea');
      else
        disp('     Spectrum for pure wind sea');
      end
    end
  else %swell dominated seas
   
    % Primary peak (swell)
    Ns  = K0*sqrt(Hm0)+K00;             %high frequency exponent
    Ms  = M0;                           %spectral width exponent
    Rps = min((1-A20)*exp(-(Eu/A2)^2)+A20,1);
    Hps = Rps*Hm0;                      % significant waveheight swell
    Tps = Tp;                           % primary peak period
    % peak enhancement factor
    gammas = KG*(1+KG0*exp(-Hm0/KG1))*(2*pi/gravity1*Hm0/(Tf^2))^r*(1+A3*Eu);
    gammas = max(gammas,1);
  
    % Secondary peak (wind)
    Nw   = Ns;                       % high frequency exponent
    Mw   = M0*(1-B2*exp(-Hm0/B3));   % spectral width exponent
    Rpw  = sqrt(1-Rps^2);
    Hpw  = Rpw*Hm0;                  % significant waveheight wind
   
    C = (Nw-1)/Mw;
    B = Nw/Mw;
    G0w    = B^C*Mw/gamma(C);%normalizing factor
    %G0w = exp(C*log(B)+log(Mw)-gammaln(C))
    %G0w  = Mw/((B)^(-C)*gamma(C)); 

    if Hpw>0,
      Tpw  = (16*S0*(1-exp(-Hm0/S1))*(0.4)^Nw/( G0w*Hpw^2))^(-1/(Nw-1));
    else
      Tpw = inf;
    end
    %Tpw  = max(Tpw,2.5)
    gammaw = 1;
    if monitor
      if Rpw > 0.1
        disp('     Spectrum for swell dominated sea');
      else
        disp('     Spectrum for pure swell sea');
      end
    end
  end
  
  if monitor,
    if (3.6*sqrt(Hm0)<= Tp && Tp<=5*sqrt(Hm0))
      disp('     Jonswap range');
    end
    disp(['Hm0 = ' num2str(Hm0)]);
    disp(['Ns, Ms = ' num2str([Ns Ms]) '  Nw, Mw = ' num2str([Nw Mw])]);
    disp(['gammas = ' num2str(gammas) ' gammaw = ' num2str(gammaw)]);
    disp(['Rps = ' num2str(Rps) ' Rpw = ' num2str(Rpw)]);
    disp(['Hps = ' num2str(Hps) ' Hpw = ' num2str(Hpw)]);
    disp(['Tps = ' num2str(Tps) ' Tpw = ' num2str(Tpw)]);
  end

  %G0s=Ms/((Ns/Ms)^(-(Ns-1)/Ms)*gamma((Ns-1)/Ms )); %normalizing factor

  % Wind part
  Hwind  = mkjonswap('Hm0',Hpw,'Tp',Tpw,'gamma',gammaw,'N',Nw,'M',Mw,'method',options.method,'chkseastate','off');
  % Swell part
  Hswell = mkjonswap('Hm0',Hps,'Tp',Tps,'gamma',gammas,'N',Ns,'M',Ms,'method',options.method,'chkseastate','off');

end % mktspec

