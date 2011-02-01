function [x2,x,svec,dvec,A]=spec2nlsdat(S,np,dt,iseed,method,truncationLimit)
%SPEC2NLSDAT Simulates a Randomized 2nd order non-linear wave X(t) 
%            given the spectral density S. 
%           
%   CALL: [xs2 xs1] = spec2nlsdat(S,[np cases],dt,iseed,method,fnLimit);
% 
%   xs2   = a cases+1 column matrix  ( t,X1(t) X2(t) ...). 
%          (1'st + 2'nd order components)
%   xs1   = a cases+1 column matrix  ( t,X1(t) X2(t) ...). 
%          (1'st order component)
%   S     = a spectral density structure
%   np    = giving np load points.  (default length(S)-1=n-1).
%           If np>n-1 it is assummed that S(k)=0 for all k>n-1
%   cases = number of cases (default=1) 
%   dt    = step in grid (default dt is defined by the Nyquist freq)
%   iseed = starting seed number for the random number generator 
%          (default none is set)
%  method = 'apStochastic'    : Random amplitude and phase (default)
%           'aDeterministic'  : Deterministic amplitude and random phase
%           'apDeterministic' : Deterministic amplitude and phase
% fnLimit = normalized upper frequency limit of spectrum for 2'nd order
%           components. The frequency is normalized with 
%           sqrt(gravity*tanh(kbar*waterDepth)/Amax)/(2*pi)
%           (default sqrt(2), i.e., Convergence criterion).
%           Other possible values are:
%            sqrt(1/2)  : No bump in trough criterion
%            sqrt(pi/7) : Wave steepness criterion   
%  
%  SPEC2NLSDAT performs a Fast simulation of Randomized 2nd order non-linear 
%  waves by summation of sinus functions with random amplitudes and 
%  phase angles.  The extent to which the simulated result are applicable
%  to real seastates are dependent on the validity of the assumptions:
%
%  1) Seastate is unidirectional
%  2) the surface elevation is adequately represented by 2nd order random
%     wave theory
%  3) The first order component of the surface elevation is a Gaussian
%     random process.
%
%  NOTE :  If the spectrum does not decay rapidly enough towards zero, the
%  contribution from the 2nd order wave components at the upper tail can
%  be very large and unphysical.
%  To ensure convergence of the perturbation series, the upper tail of the
%  spectrum is truncated at FNLIMIT in the calculation of the 2nd order
%  wave components, i.e., in the calculation of sum and difference
%  frequency effects. This may also be combined with the elimination of
%  second order effects from the spectrum, i.e., extract the linear
%  components from the spectrum. One way to do this is to use SPEC2LINSPEC.   
%
% Example:
%  np =100; dt = .2;
%  [x1, x2] = spec2nlsdat(jonswap,np,dt);
%  waveplot(x1,'r',x2,'g',1,1)  
% 
%  %More extensive test
%  Sj = jonswap
%  [x2, x1] = spec2nlsdat(Sj,[20000,20]);
%  [sk,ku]= spec2skew(Sj);
%  truth1 = [0, sqrt(spec2mom(Sj,1)), sk, ku]; 
%  funs = {@mean,  @std, @skew, @kurt}
%  for i = 1:4,
%      trueval = truth1(i);
%      fun = funs{i};
%      res = fun(x2(:,2:end), 1);
%      m = mean(res);
%      sa = std(res);
%      [abs(m-trueval)<sa, trueval, m, sa]
%   end
%
% See also  spec2linspec, spec2sdat, cov2sdat

% Reference 
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
% Revised pab Feb2004
% - changed seed to state   
% revised pab 25Jan2004
%  - changed the truncation at the upper tail of the spectrum  
%  - added truncationLimit to input
% revised pab 11Nov2003
%  changed call from disufq1 to disufq  
% revised pab 22.07.2002
% revised pab 15.03.2002
% -new call to disufq
% - added nargchk
% by pab 21.01.2001

% TODO % Check the methods: 'apdeterministic' and 'adeterministic' 
  
% Variables controlling the truncation of the spectrum for sum and
% difference frequency effects  
reltol2ndorder     = 1e-3; %
%truncationLimit = 1.5;

error(nargchk(1,6,nargin))

ftype = freqtype(S); %options are 'f' and 'w' and 'k'
n     = length(S.(ftype));

numWaves = 1000; % Typical number of waves in 3 hour seastate
%C = pi/7; % Wave steepness criterion
%C = 1/2 ; % No bump in trough
C = 2;    % Convergence criterion as given in Nestegaard and Stokka (1995)
if (nargin<6 || isempty(truncationLimit)), 
  truncationLimit = sqrt(C);
end
if (nargin<5 || isempty(method)), method = 'apstochastic'; end
if (nargin>3 && ~isempty(iseed)),
  try
    randn('state',iseed);
  catch
    randn('seed',iseed);
  end
end  % set the the seed  
if (nargin<2 || isempty(np)),     np = n-1; end
if (nargin>2 && ~isempty(dt)),    S = specinterp(S,dt);end  % interpolate spectrum  
							  
switch  length(np) 
  case 1, cases=1; 
  case 2, cases=np(2); np=np(1);
  otherwise, error('Wrong input. Too many arguments')
end
np = np + mod(np,2); % make sure np is even    

fs    = S.(ftype);
Si    = S.S(2:end-1);
h     = S.h;
if isempty(h), h = inf;end


switch ftype
  case 'f'   
  case {'w','k'}
    Si = Si*2*pi;
    fs = fs/2/pi;
  otherwise
    error('Not implemented for wavenumber spectra')
end
dT = 1/(2*fs(end)); % dT



df = 1/(np*dT);


% interpolate for freq.  [1:(N/2)-1]*df and create 2-sided, uncentered spectra
% ----------------------------------------------------------------------------
f = (1:(np/2)-1)'*df;

fs(1)   = []; 
fs(end) = [];
Fs      = [0; fs(:); (np/2)*df];
Su      = [0; abs(Si(:))/2; 0];

Smax = max(Su);
%Si = interp1(Fs,Su,f,'linear');
Si = interp1q(Fs,Su,f);

% If the spectrum does not decay rapidly enough towards zero, the
% contribution from the wave components at the  upper tail can be very
% large and unphysical.
% To ensure convergence of the perturbation series, the upper tail of the
% spectrum is truncated in the calculation of sum and difference
% frequency effects.
% Find the critical wave frequency to ensure convergence. 
tmp  = find(Si>Smax*reltol2ndorder);
switch 2
 case 1
  nmax = max(tmp)+1;
  nmin = min(tmp)+1;
 case 2,
  Hm0  = spec2char(S,'Hm0');
  Tm02 = spec2char(S,'Tm02');
  waterDepth = abs(S.h);
  kbar = w2k(2*pi/Tm02,0,waterDepth);
  
  Amax = sqrt(2*log(numWaves))*Hm0/4; % Expected maximum amplitude for 1000 waves seastate
  
  fLimitUp = truncationLimit*sqrt(gravity*tanh(kbar*waterDepth)/Amax)/(2*pi);
  fLimitLo = sqrt(gravity*tanh(kbar*waterDepth)*Amax/waterDepth^3)/(2*pi);
  
  nmax   = min(max(find(f<=fLimitUp)),max(find(Si>0)))+1;
  nmin   = max(min(find(fLimitLo<=f)),min(tmp))+1;
  %nmin = min(tmp)+1;
end
if isempty(nmax),nmax = np/2;end
if isempty(nmin),nmin = 2;end % Must always be greater than 1
fLimitUp = df*nmax;
fLimitLo = df*nmin;

disp(sprintf('2nd order frequency Limits = %g,%g',fLimitLo, fLimitUp))

Su = [0; Si; 0; Si((np/2)-1:-1:1)];

clear Si Fs

T       = (np-1)*dT;
x       = zeros(np,cases+1);
x(:,1)  = linspace(0,T,np)'; %(0:dT:(np-1)*dT).';
x2      = x;

w  = 2*pi*[0; f; np/2*df];
g  = gravity;
kw = w2k(w ,[],h,g);

% Generate standard normal random numbers for the simulations
% -----------------------------------------------------------
Zr = randn((np/2)+1,cases);
Zi = [zeros(1,cases); randn((np/2)-1,cases); zeros(1,cases)];

A                = zeros(np,cases);
A(1:(np/2+1),:)  = Zr - sqrt(-1)*Zi; clear Zr Zi
A((np/2+2):np,:) = conj(A(np/2:-1:2,:));
A(1,:)           = A(1,:)*sqrt(2);
A((np/2)+1,:)    = A((np/2)+1,:)*sqrt(2);

% Make simulated time series
% --------------------------

%mean(abs(A(:)))
Ssqr = sqrt(Su*df/2);
%max(Ssqr)
if strncmpi(method,'apdeterministic',3)
  % Deterministic amplitude and phase
  A(2:(np/2),:)    = A(2,1);
  A((np/2+2):np,:) = conj(A(2,1)); 
  A = sqrt(2)*Ssqr(:,ones(1,cases)).*exp(sqrt(-1)*atan2(imag(A),real(A)));
elseif strncmpi(method,'adeterministic',3)
   % Deterministic amplitude and random phase
  A = sqrt(2)*Ssqr(:,ones(1,cases)).*...
      exp(sqrt(-1)*atan2(imag(A),real(A)));
else
   % stochastic amplitude and phase
  A = A.*Ssqr(:,ones(1,cases));
end
%max(abs(A))
clear Su Ssqr
  

x(:,2:end) = real(fft(A));
  
if nargout>3, 
   %compute the sum and frequency effects separately
  [svec, dvec] = disufq((A.'),w,kw,min(h,10^30),g,nmin,nmax);
  svec = svec.';
  dvec = dvec.';
  
  x2s  = fft(svec); % 2'nd order sum frequency component 
  x2d  = fft(dvec); % 2'nd order difference frequency component
  
  % 1'st order + 2'nd order component.
  x2(:,2:end) =x(:,2:end)+ real(x2s(1:np,:))+real(x2d(1:np,:)); 
else
  svec = disufq((A.'),w,kw,min(h,10^30),g,nmin,nmax).';
  
  x2o  = fft(svec); % 2'nd order component 
  
  
  % 1'st order + 2'nd order component.
  x2(:,2:end)=x(:,2:end)+ real(x2o(1:np,:)); 
end
return




