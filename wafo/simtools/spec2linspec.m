function [SL,SN]=spec2linspec(S,np,dt,iseed,fnLimit)
%SPEC2LINSPEC  Separates the linear component of the Spectrum 
%              according to 2nd order wave theory
%           
%   CALL: [SL,SN] = spec2linspec(S,[np cases],dt,iseed,fnLimit);
% 
%   SL    = spectral density structure with linear components only.
%   SN    = spectral density structure with non-linear components only.
%   S     = spectral density structure with linear and non-linear.
%           components
%   np    = giving np load points.  (default length(S)-1=n-1).
%           If np>n-1 it is assummed that S(k)=0 for all k>n-1
%   cases = number of cases (default=20) 
%   dt    = step in grid (default dt is defined by the Nyquist freq)
%   iseed = starting seed number for the random number generator 
%          (default none is set)
% fnLimit = normalized upper frequency limit of spectrum for 2'nd order
%           components. The frequency is normalized with 
%           sqrt(gravity*tanh(kbar*waterDepth)/Amax)/(2*pi)
%           (default sqrt(2), i.e., Convergence criterion).
%           Generally this should be the same as used in the final
%           non-linear simulation (see example below).
%
%  SPEC2LINSPEC separates the linear and non-linear component of the Spectrum 
%  according to 2nd order wave theory. This is useful when simulating
%  non-linear waves because:
%  If the spectrum does not decay rapidly enough towards zero, the
%  contribution from the 2nd order wave components at the upper tail can
%  be very large and unphysical.
%  Another option to ensure convergence of the perturbation series in the
%  simulation, is to truncate the upper tail of the
%  spectrum at FNLIMIT in the calculation of the 2nd order
%  wave components, i.e., in the calculation of sum and difference
%  frequency effects. 
%
% Example:
%  np = 10000;
%  iseed = 1;
%  pflag = 2;
%  S  = jonswap(10);
%  fnLimit = inf;  
%  [SL,SN] = spec2linspec(S,np,[],[],fnLimit);
%  x0 = spec2nlsdat(SL,8*np,[],iseed,[],fnLimit);
%  x1 = spec2nlsdat(S,8*np,[],iseed,[],fnLimit); 
%  x2 = spec2nlsdat(S,8*np,[],iseed,[],sqrt(2));  
%  Se0 = dat2spec(x0);
%  Se1 = dat2spec(x1);
%  Se2 = dat2spec(x2); 
%  clf  
%  plotspec(SL,'r',pflag),  % Linear components
%  hold on
%  plotspec(S,'b',pflag)    % target spectrum for simulated data
%  plotspec(Se0,'m',pflag), % approx. same as S 
%  plotspec(Se1,'g',pflag)  % unphysical spectrum
%  plotspec(Se2,'k',pflag)  % approx. same as S
%  axis([0 10 -80 0])
%  hold off
%  
% See also  spec2nlsdat

% Reference 
% P. A. Brodtkorb (2004), 
% The probability of Occurrence of dangerous Wave Situations at Sea.
%  Dr.Ing thesis, Norwegian University of Science and Technolgy, NTNU,
%  Trondheim, Norway.
%  
% Nestegaard, A  and Stokka T (1995)
% A Third Order Random Wave model.
% In proc.ISOPE conf., Vol III, pp 136-142.
%
% R. S Langley (1987)
% A statistical analysis of non-linear random waves.
% Ocean Engng, Vol 14, pp 389-407
%
% Marthinsen, T. and Winterstein, S.R (1992)
% 'On the skewness of random surface waves'
% In proc. ISOPE Conf., San Francisco, 14-19 june.

% tested on: Matlab 5.3
% History:
% Revised pab Nov 2004  
% changed the default constant controlling its
%                 performance. Can be improved further
% by pab 13.08.2002

% TODO % Replace inputs with options structure
% TODO % Can be improved further.
error(nargchk(1,5,nargin))

% Define some constants
%fnLimit = sqrt(inf)
method = 'apstochastic';
trace     = 1; % trace the convergence
maxSim    = 30;
tolerance = 5e-4;
cases     = 20;
L         = 200; %maximum lag size of the window function used in
                 %spectral estimate
ftype = freqtype(S); %options are 'f' and 'w' and 'k'
n     = length(S.(ftype));
switch ftype
 case 'f', 
  ftype = 'w';
  S = ttspec(S,ftype);
end
Hm0  = spec2char(S,'Hm0');
Tm02 = spec2char(S,'Tm02');

if (nargin<5 || isempty(fnLimit))
  fnLimit = sqrt(2);
end
if (nargin>3 && ~isempty(iseed)), 
  randn('state',iseed); % set the the seed 
else
  iseed = 0;
end   
if (nargin<2 || isempty(np)),     np = max(n-1,5000);  end
if (nargin>2 && ~isempty(dt)),    S = specinterp(S,dt);end  % interpolate spectrum  
							  
switch  length(np) 
  case 1, 
  case 2, cases=np(2); np=np(1);
  otherwise, error('Wrong input. Too many arguments')
end
np = np + mod(np,2); % make sure np is even    



waterDepth = abs(S.h);
kbar = w2k(2*pi/Tm02,0,waterDepth);

% Expected maximum amplitude for 1000 waves seastate
numWaves = 10000;  
Amax = sqrt(2*log(numWaves))*Hm0/4; 
  
fLimitLo = sqrt(gravity*tanh(kbar*waterDepth)*Amax/waterDepth^3);


freq = S.(ftype);
freq(end) = freq(end)-sqrt(eps);
Hw2  = 0;

SL = S;

indZero = find(freq<fLimitLo);
if any(indZero)
  SL.S(indZero) = 0;
end
maxS = max(S.S);
%Fs = 2*freq(end)+eps; % sampling frequency

for ix=1:maxSim
  [x2,x1] = spec2nlsdat(SL,[np,cases],[],iseed,method,fnLimit);
  %x2(:,2:end) = x2(:,2:end) -x1(:,2:end);
  S2 = dat2spec(x2,L);   
  S1 = dat2spec(x1,L);
  %[tf21,fi] = tfe(x2(:,2),x1(:,2),1024,Fs,[],512);
  %Hw11 = interp1q(fi,tf21.*conj(tf21),freq);
  if 1 || ix==1
    Hw1  = exp(interp1( S2.w,log(abs(S1.S./S2.S)),freq, 'linear'));  
  else
    % Geometric mean
    Hw1 =  exp((interp1( S2.w,log(abs(S1.S./S2.S)),freq, 'linear')+log(Hw2))/2);  
  end
  %Hw1  = (interp1q( S2.w,abs(S1.S./S2.S),freq)+Hw2)/2;
  %plot(freq, abs(Hw11-Hw1),'g')
  %title('diff')
  %pause
  %clf
  
  %d1 = interp1q( S2.w,S2.S,freq);;
  
  SL.S = (Hw1.*S.S);
  
  if any(indZero)
    SL.S(indZero) = 0;
  end
  k = find(SL.S< 0);
  if any(k), % Make sure that the current guess is larger than zero
    %k
    %Hw1(k)
    Hw1(k)  = min((S1.S(k)*0.9),(S.S(k)));
    SL.S(k) = max(Hw1(k).*S.S(k),eps);
  end
  Hw12 = Hw1-Hw2;
  maxHw12 = max(abs(Hw12));
  if trace==1,
    figure(1), plot(freq,Hw1,'r'), hold on, title('Hw')
    set(gca,'yscale','log')
    figure(2), plot(freq,abs(Hw12),'r'), hold on, title('Hw-HwOld')
    set(gca,'yscale','log')
    drawnow
    pause(3)
    figure(1), plot(freq,Hw1,'b'), hold on, title('Hw')
    figure(2), plot(freq,abs(Hw12),'b'), hold on, title('Hw-HwOld')
    figtile
  end
   
   disp(['Iteration ',num2str(ix),...
        ' Hw12 :  ' num2str(maxHw12), ...
	' Hw12/maxS : ' num2str(maxHw12/maxS)]),
  if (maxHw12<maxS*tolerance  && Hw1(end)<Hw2(end) )
    break
  end
  Hw2 = Hw1;
end

%Hw1(end)
%maxS*1e-3
%if Hw1(end)*S.>maxS*1e-3,
%  warning('The Nyquist frequency of the spectrum may be too low')
%end

SL.date = datestr(now);
if nargout>1
  SN   = SL;
  SN.S = S.S-SL.S;
  SN.note = cat(2,SN.note,' non-linear component (spec2linspec)');
end
SL.note = cat(2,SL.note,' linear component (spec2linspec)');

return



