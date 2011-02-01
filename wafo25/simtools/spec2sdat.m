function [x,xder]=spec2sdat(S,np,dt,iseed,method)
%SPEC2SDAT Simulates a Gaussian process and its derivative from spectrum
%           
%   CALL: [xs, xsder] = spec2sdat(S,[np cases],dt,iseed,method);
% 
%         xs    = a cases+1 column matrix  ( t,X1(t) X2(t) ...).
%         xsder = a cases+1 column matrix  ( t,X1'(t) X2'(t) ...). 
%         S     = a spectral density structure
%         np    = giving np load points.  (default length(S)-1=n-1).
%                 If np>n-1 it is assummed that R(k)=0 for all k>n-1
%         cases = number of cases, i.e. number of replicates (default=1) 
%            dt = step in grid (default dt is defined by the Nyquist freq)
%         iseed = starting seed number for the random number generator 
%                 (default none is set)
%        method = 'exact'  : simulation using cov2sdat 
%                 'random' : random phase and amplitude simulation (default)
%
%  SPEC2SDAT performs a fast and exact simulation of stationary zero mean 
%  Gaussian process through circulant embedding of the covariance matrix
%  or by summation of sinus functions with random amplitudes and random
%  phase angle.
%  
%  If the spectrum has a non-empty field .tr, then the transformation is 
%  applied to the simulated data, the result is a simulation of a transformed
%  Gaussian process.
%
%  NB! The method 'exact' simulation may give high frequency ripple when 
%  used with a small dt. In this case the method 'random' works better. 
%
% Example:
%  np =100; dt = .2;
% [x1 x2] = spec2sdat(jonswap,np,dt);
% waveplot(x1,'r',x2,'g',1,1)  
%
%  %More extensive test
%  Sj = jonswap
%  [x2, x1] = spec2sdat(Sj,[20000,20]);
%  [sk,ku]= spec2skew(Sj);
%  truth1 = [0, sqrt(spec2mom(Sj,1)), sk, ku-3]; 
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
% See also  cov2sdat, gaus2dat

% Reference 
% C.R Dietrich and G. N. Newsam (1997)
% "Fast and exact simulation of stationary 
%  Gaussian process through circulant embedding 
%  of the Covariance matrix" 
%  SIAM J. SCI. COMPT. Vol 18, No 4, pp. 1088-1107
%
% Hudspeth, R.T. and Borgman, L.E. (1979)
% "Efficient FFT simulation of Digital Time sequences"
% Journal of the Engineering Mechanics Division, ASCE, Vol. 105, No. EM2, 
% 

% Tested on: Matlab 5.3
% History:
% Revised pab Mar2005
% 
% Revised pab Feb2004
% - changed seed to state   
% revised pab 21.07.2002
% - 
% revised pab 12.10.2001
% - added example
% - the derivative, xder, is now calculated 
%   correctly using fft for method = 'random'.
%   derivate.m is no longer needed!
% revised pab 11.10.2001
% - added a trick to avoid adding high frequency noise to spectrum when
%   method = 'exact' 
% revised pab 21.01.2001
% - random is now available for 1D wavenumber spectra as well
% revised ir 11.06.00 - introducing dt
% revised es 24.05.00 - fixed default value for np, and small changes to help
% revised pab 13.03.2000
% - fixed default value for np
% revised pab 24.01.2000
% - added method random from L. Borgman fwavsim (slightly faster and needs 
%   less memory)
% - added iseed
% revised es 12.01.2000:  enable dir. spectrum input (change of spec2cov call)
% revised pab 12.10.1999
%  simplified call by calling cov2sdat
%    last modified by Per A. Brodtkorb 19.08.98

if nargin<5||isempty(method)
  method='random';% 'exact' is slightly slower than method=='random'
end
if nargin<4||isempty(iseed)
 iseed=[];
else
  try
    randn('state',iseed)
  catch
    randn('seed',iseed); % set the the seed  
  end
end

if nargin<3||isempty(dt)
 S1=S;
else
  S1=specinterp(S,dt); % interpolate spectrum  
end

ftype = freqtype(S);
freq  = S.(ftype);

Nt = length(freq);

S = S1;

if nargin<2||isempty(np)
  np=Nt-1;
end
if strcmpi(method,'exact') 
  R = spec2cov(S,0,[],1,0,0);
  
  if strcmpi(ftype,'f')
    T = Nt/(2*freq(end));
  else
    T = Nt*pi/freq(end);
  end
  ltype = lagtype(R);

  ix = find(R.(ltype)>T,'first');
  
  % Trick to avoid adding high frequency noise to the spectrum
  if ~isempty(ix),
    R.R(ix:end)=0; 
  end
  
  if nargout>1
    [x, xder]=cov2sdat(R,np,iseed);
  else
    x=cov2sdat(R,np,iseed);
  end
  return
end
   
switch  length(np) 
  case 1, cases=1; 
  case 2, cases=np(2); np=np(1);
  otherwise, error('Wrong input. Too many arguments')
end    
if mod(np,2),
  np=np+1;
end % make sure it is even


if 0
  % Do the simulation with normalized spectrum
  % Normalize the spectrum
  [Sn,mn4,m0,m2] = specnorm(S);
  % Normalization constants
  tnorm = 2*pi*sqrt(m0/m2);
  xnorm = sqrt( tnorm*m0/(2*pi));
  fs    = Sn.(ftype);
  Si    = Sn.S(2:end-1);
else
  tnorm = 1;
  xnorm = 1;
  
  fs    = S.(ftype);
  Si    = S.S(2:end-1);
  
  switch ftype
    case 'f'   
    case {'w','k'}
      Si=Si*2*pi;
      fs=fs/2/pi;
    otherwise
      error('Not implemented for wavenumber spectra')
  end
end

x=zeros(np,cases+1);
if nargout==2
  xder=x;
end


dT = 1/(2*fs(end)); % dT

df = 1/(np*dT);


% interpolate for freq.  [1:(N/2)-1]*df and create 2-sided, uncentered spectra
% ----------------------------------------------------------------------------
f = (1:(np/2)-1).'*df;

fs(1)   = []; 
fs(end) = [];
Fs      = [0; fs(:); (np/2)*df];
Su      = [0; abs(Si(:))/2; 0];

%Si = interp1(Fs,Su,f,'linear');
Si = interp1q(Fs,Su,f);
Su=[0; Si; 0; Si((np/2)-1:-1:1)];

clear Si Fs

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

T    = (np-1)*dT;
Ssqr = sqrt(Su*df/2);

% stochastic amplitude
A    = A.*Ssqr(:,ones(1,cases));


% Deterministic amplitude 
%A     =  sqrt(2)*Ssqr(:,ones(1,cases)).*exp(sqrt(-1)*atan2(imag(A),real(A)));
clear Su Ssqr
  

x(:,2:end) = xnorm*real(fft(A));

x(:,1)     = linspace(0,T*tnorm,np)'; %(0:dT:(np-1)*dT).';




if nargout==2, 
  %xder=derivate(x(:,1),x(:,2:end));
  
  % new call pab 12.10.2001
  w = 2*pi*[0; f;0;-f(end:-1:1)] ;
  A = -i*A.*w(:,ones(1,cases));  
  xder(:,2:(cases+1)) = real(fft(A))*xnorm/tnorm;
  xder(:,1)           = x(:,1);
end


if isfield(S,'tr') && ~isempty(S.tr)
  disp('   Transforming data.')
  g=S.tr;
  G=fliplr(g); % the invers of g
  if nargout==2, % gaus2dat
    for ix=1:cases
      tmp=tranproc([x(:,ix+1) xder(:,ix+1)],G); 
      x(:,ix+1)=tmp(:,1);
      xder(:,ix+1)=tmp(:,2);
    end
  else
    for ix=1:cases % gaus2dat
      x(:,ix+1)=tranproc(x(:,ix+1),G); 
    end
  end
end



return

% Old call kept just in case



% df=1/(np*dT);
% 
% 
% % Interpolate for freq.  [1:(N/2)-1]*df and create 2-sided, uncentered spectra
% % ----------------------------------------------------------------------------
% f=(1:(np/2)-1)'*df;
% 
% fs(1)=[];fs(end)=[];
% Fs=[0; fs(:); (np/2)*df];
% Su=[0; Si(:); 0];
% Si = interp1(Fs,Su,f,'linear')/2;
% 
% Su=[0; Si; 0; Si((np/2)-1:-1:1)];
% clear Si
% 
% % Generate standard normal random numbers for the simulations
% % -----------------------------------------------------------
% %Zr=randn((np/2)+1,cases);
% %Zi=[zeros(1,cases); randn((np/2)-1,cases); zeros(1,cases)];
% %AA=Zr-sqrt(-1)*Zi;
% 
% AA = randn((np/2)+1,cases)-sqrt(-1)*[zeros(1,cases); randn((np/2)-1,cases); zeros(1,cases)];
% A=[AA;  conj(AA(np/2:-1:2,:))];
% clear AA 
% A(1,:)=A(1,:)*sqrt(2);
% A((np/2)+1,:)=A((np/2)+1,:)*sqrt(2);
% 
% % Make simulated time series
% % --------------------------
% T    = np*dT;
% Ssqr = sqrt(Su*df/2);
% A    = A.*Ssqr(:,ones(1,cases));
% clear Su Ssqr
% 
% x(:,2:end)=real(fft(A));
% 
% x(:,1)=(0:dT:(np-1)*dT).';