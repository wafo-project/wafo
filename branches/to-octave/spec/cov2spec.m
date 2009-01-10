function S = cov2spec(R,rate,nugget,trunc,fast)
%COV2SPEC Computes spectral density given the auto covariance function 
%
%  CALL: S = cov2spec(R,rate,nugget,trunc,fast);
%
%       R    = a covariance structure 
%
%       S    = a spectral density structure 
%
%  Optional parameters:
%
%       rate = 1,2,4,8...2^r, interpolation rate for f
%              (default = 1, no interpolation)
%
%     nugget = add a nugget effect to ensure that round off errors
%              do not result in negative spectral estimates
%              (default 0) Good choice might be 10^-12.
%
%     trunc  = truncates all spectral values where S/max(S) < trunc
%              0 <= trunc <1   This is to ensure that high frequency 
%              noise is not added to the spectrum.  (default 1e-5)
%
%     fast   = 1 if zero-pad to obtain power of 2 length R.R (default)
%              0 if no zero-padding of R.R, slower but more accurate.  
%
% NB! This routine requires that the covariance is evenly spaced
%     starting from zero lag. Currently only capable of 1D matrices.
%
% Example:
%  R = createcov;
%  L = 129;  
%  R.t = linspace(0,75,L).';
%  R.R = zeros(L,1);
%  win = parzen(40);  
%  R.R(1:20) = win(21:40);
%  S0 = cov2spec(R);
%  S = jonswap;
%  S1 = cov2spec(spec2cov(S));
% assert(all(abs(S1.S-S.S)<1e-4) ,'COV2SPEC')
%  
% See also  spec2cov, datastructures
  
% tested on: matlab 5.3
% history:
% Revised pab   
% -Increased the precision in calculations
% - deleted L from input, added fast instead
% revised jr 11.07.2000
%  - line 158. Changed  n  to  nf .
% revised pab 13.03.2000
%  - commented out warning messages
% revised pab 16.01.2000
%  - fixed a bug: truncate negative values of the spectral density to
%     zero to make sure these spurious points are not added to the
%     spectrum
%  - added trunc: this is a fix: truncating the values of the spectral density to
%     zero if S/max(S) < trunc. This is to ensure that high frequency 
%     noise is not added to the spectrum.
% revised by pab 23.09.1999  
  
%        dT   = time-step between data points.(default = T(2)-T1)).

% add a nugget effect to ensure that round off errors
% do not result in negative spectral estimates

if ~isstruct(R)
  error('Incorrect input Covariance, see help datastructures')
end
if ndims(R.R)>2||(min(size(R.R))>1),
  error('This function is only capable of 1D covariances')
end

if nargin<3||isempty(nugget)
  nugget = 0;%10^-12;
end

if nargin<4||isempty(trunc)
  trunc = 1e-5;
else
  trunc = min(abs(trunc),1); % make sure it
end
if nargin<5 || isempty(fast)
  fast = 1;
end

%comment=0;



n     = length(R.R);
names = fieldnames(R);
ind   = find(strcmpi(names,'x')+strcmp(names,'y')+strcmp(names,'t')); 
        %options are 'x' and 'y' and 't'

lagtype = lower(names{ind});
if strcmpi(lagtype,'t')
  spectype = 'freq';
  ftype     = 'w'; %options f=S(f), w=S(w)
else
  spectype = 'k1d';
  ftype     = 'k';  %options f=S(f), w=S(w)
end
ti = R.(lagtype);
dT = ti(2)-ti(1); % ti 



if nargin<2||isempty(rate)
  rate = 1;%interpolation rate
else
  rate = 2^nextpow2(rate);%make sure rate is a power of 2
end

% add a nugget effect to ensure that round off errors
% do not result in negative spectral estimates
ACF    = R.R(:);
ACF(1) = R.R(1) +nugget;
% embedding a circulant vector and Fourier transform
if fast
  nfft = 2^nextpow2(2*n-2);
else
  nfft = 2*n-2;
end
nf   = nfft/2;% number of frequencies

ACF  = [ACF;zeros(nfft-2*n+2,1);ACF(n-1:-1:2)];

Rper = real(fft(ACF,nfft));% periodogram

k = find(Rper<0);
if any(k)
 % disp('Warning: negative spectral estimates')
 % disp('Apply the parzen windowfunction ') 
 % disp('to the ACF in order to avoid this.')
 % disp('The returned result is now only an approximation.')
  %min(Rper(k))
  Rper(k)=0; % truncating negative values to ensure that 
            % that high frequency noise is not added to 
	    % the Spectrum
end


k = find(Rper/max(Rper(:))<trunc);
if any(k)
  Rper(k)=0; % truncating small values to ensure that 
            % that high frequency noise is not added to 
	    % the Spectrum
end

%size(nf),size(dT)
S      = createspec(spectype,ftype);
S.tr   = R.tr;
S.h    = R.h;
S.norm = R.norm;
S.note = R.note;
%fn = linspace(0,1,
switch ftype
case {'w'}
   S.w = (0:(nf))'/nf*pi/dT;  % (rad/s)
  S.S  = abs(Rper(1:(nf+1)))*dT/pi; % (m^2*s/rad) one-sided spectrum
case {'f'} % spectype == sf
  S.f =(0:(nf))'/nf/2/dT;      % frequency Hz if dT is in seconds
  S.S =2*abs(Rper(1:(nf+1)))*dT; % (m^2*s) one sided spectrum
case {'k'}
  S.k =(0:(nf))'/nf*pi/dT;% (rad/m) 
  S.S = abs(Rper(1:(nf+1)))*dT/pi; % (m^3/rad) one-sided wavenumber spectrum
end

if rate>1
  w = S.(ftype);
  wi = linspace(w(1),w(end),nf*rate).';
  S.S = interp1(w,S.S,wi);
  S.(ftype) = wi;
end





