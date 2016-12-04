function [S, Sw, Ss]=torsethaugen(w1,sdata,plotflag)
% TORSETHAUGEN Calculates a double peaked (swell + wind) spectrum 
%
%  CALL: [S, Sw, Ss] = torsethaugen(w,data,plotflag); 
%        [S, Sw, Ss] = torsethaugen(wc,data,plotflag); 
%
%   S, Ss, Sw = spectrum struct for total, swell and wind, respectively.  
%        w    = angular frequency (default linspace(0,wc,257))
%        wc   = angular cutoff frequency (default 33/Tp)
%        data = [Hm0 Tp A]
%               Hm0 = Significant wave height      (default 7)
%               Tp  = 2*pi/wp, primary peak period (deafult 11)
%               A   = alpha, normalization factor, (default -1)
%                   A<0  : A calculated by integration so that int S dw =Hm0^2/16
%                   A==0 : A = (1+f1(N,M)*log(gammai)^f2(N,M))/gammai, original 
%                         parametric normalization  
%    plotflag = 0, do not plot the spectrum (default).
%               1, plot the spectrum.
%
% For zero values, NaN's or parameters not specified in DATA the
% default values are used. 
%
% Model:  S(w)=Ss(w)+Sw(w) 
% where Ss and Sw are modified JONSWAP spectrums for 
% swell and wind peak, respectively.
% The energy is divided between the two peaks according
% to empirical parameters, which peak that is primary depends on parameters.   
% The empirical parameters are found for classes of Hm0 and Tp,
% originating from a dataset consisting of 20 000 spectra divided
% into 146 different classes of Hm0 and Tp. (Data measured at the
% Statfjord field in the North Sea in a period from 1980 to 1989.)
% The range of the measured  Hm0 and Tp for the dataset
% are from 0.5 to 11 meters and from 3.5 to 19 sec, respectively.
% See Torsethaugen (1996).
%
% Preliminary comparisons with spectra from other areas indicate that 
% some of the empirical parameters are dependent on geographical location.
% Thus the model must be used with care for other areas than the
% North Sea and sea states outside the area where measured data 
% are available.
%
% Example:
% [S,Sw,Ss]=torsethaugen([],[6 8],1);
%
% close all
%
% See also  jonswap, pmspec.

% References 
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


% Tested on Matlab 5.3 
% History: 
% revised pab april 2007
% -replaced code with a call to mktorsethaugen
% Revised pab april 2005
% revised pab 12.12.2000
% - fixed a bug: Hpw = 0 or Hps = 0 is now handled properly
% - sdata may now contain NaN's 
% revised pab 16.02.2000
%  -added ih, see at the end
% revised pab 20.12.1999
%  -added norm
%  -changed default value of Hm0 to 7 (equal to jonswap)
% revised pab 01.12.1999
%  - updated the reference
% revised  es 23.09.1999
%  - updated documentation
%  - added note
% by pab 23.06.1999


% default values
Hm0   = 7;
Tp    = 11; 
Aj    = -1;
data2 = [Hm0, Tp, Aj];
nd2   = length(data2);
if nargin<3||isempty(plotflag),  plotflag = 0;       end
if (nargin>1) && ~isempty(sdata), 
  nd = length(sdata); 
  ind = find(~isnan(sdata(1:min(nd,nd2))));
  if any(ind), % replace default values with those from input data
    data2(ind) = sdata(ind);
  end
end
if (nd2>0) && (data2(1)>0), Hm0 = data2(1);end
if (nd2>1) && (data2(2)>0),  Tp = data2(2);end
if (nd2>2) && (data2(3)==0),  Aj = data2(3);end

w = [];
if nargin<1||isempty(w1), 
   wc = 33/Tp;
elseif length(w1)==1,
   wc = w1; 
else
   w = w1 ;
end
nw = 257;

if isempty(w),      w = linspace(0,wc,nw).'; end

w=w(:);

options = mktorsethaugen('defaults');
options.Hm0 = Hm0;
options.Tp = Tp;
if Aj==0
  options.method = 'parametric';
end
[tspec tspecW,tspecS] = mktorsethaugen(options); 


%n      = length(w);
S      = createspec('freq');
S.w    = w;
S.norm = 0; % not normalized spectra

Sw = S;
Sw.S = tspecW(w); % Wind part

Ss = S;
Ss.S = tspecS(w); % Swell part

S.S = Ss.S+Sw.S; % total spectrum

S.note=strcat('Torsethaugen, Hm0=',num2str(Hm0),', Tp=',num2str(Tp));

if max(Ss.S)<max(Sw.S)
   Ss.note=strcat('Swell, secondary peak of Torsethaugen, Hm0=',num2str(Hm0),', Tp=',num2str(Tp));
   Sw.note=strcat('Wind, primary peak of Torsethaugen, Hm0=',num2str(Hm0),', Tp=',num2str(Tp));
else
   Ss.note=strcat('Swell, primary peak of Torsethaugen, Hm0=',num2str(Hm0),', Tp=',num2str(Tp));
   Sw.note=strcat('Wind, secondary peak of Torsethaugen, Hm0=',num2str(Hm0),', Tp=',num2str(Tp));
end
if plotflag
  plotspec(S,plotflag,'b');
  ih = ishold;
  hold on;
  plotspec(Ss,plotflag,'b--');
  plotspec(Sw,plotflag,'b--');
  if ~ih, hold off; end
end





