function  R = spec2cov(S,nr,Nt,rate,Nx,Ny,dx,dy)
%SPEC2COV Computes covariance function and its derivatives 
%
%  CALL:  R = spec2cov(S,nr,Nt,rate,Nx,Ny,dx,dy);
%
%       R    = a covariance structure (See datastructures)
%       S    = a spectral density structure (See datastructures)
%       nr   = number of derivatives in output, nr<=4 (default = 0).
%       Nt   = number in time grid, i.e., number of time-lags
%              (default rate*(length(S.S)-1)).
%       rate = 1,2,4,8...2^r, interpolation rate for R
%               (default = 1, no interpolation) 
%     Nx,Ny   = number in space grid (default = )
%     dx,dy   = space grid step (default: depending on S)
%
% The input 'rate' gives together with the spectrum
% the t-grid-spacing: dt=pi/(S.w(end)*rate), S.w(end) is the Nyquist freq.
% This results in the t-grid: 0:dt:Nt*dt.
%
% What output is achieved with different S and choices of Nt,Nx and Ny: 
% 1) S.type='freq' or 'dir', Nt set, Nx,Ny not set: then result R(t) (one-dim)
% 2) S.type='k1d' or 'k2d', Nt set, Nx,Ny not set: then result R(x) (one-dim)
% 3) Any type, Nt and Nx set =>R(x,t); Nt and Ny set =>R(y,t)
% 4) Any type, Nt, Nx and Ny set => R(x,y,t)
% 5) Any type, Nt not set, Nx and/or Ny set => Nt set to default, goto 3) or 4)
%
% NB! This routine requires that the spectrum grid is equidistant
%     starting from zero frequency.
% NB! If you are using a model spectrum, S, with sharp edges 
%     to calculate covariances then you should probably round off the sharp
%     edges like this:
% Example:    
%    S   = jonswap; 
%    S.S([1:40 100:end]) = 0;   
%    Nt  = length(S.S)-1;  
%    R   = spec2cov(S,0,Nt);
%    win = parzen(2*Nt+1);
%    R.R = R.R.*win(Nt+1:end);
%    S1  = cov2spec(R);
%    R2  = spec2cov(S1);
%    figure(1)
%    plotspec(S),hold on, plotspec(S1,'r')
%    figure(2)
%    covplot(R), hold on, covplot(R2,[],[],'r')
%    figure(3)
%    semilogy(abs(R2.R-R.R)), hold on,
%    semilogy(abs(S1.S-S.S)+1e-7,'r')  
%  
%    close all
%
% See also  cov2spec, datastructures

% NB! requires simpson

% tested on: matlab 5.3
% history:
% -revised pab aug2007
% -updated example in help header
% -revised pab feb2007
% -replaced all call to setfield with R.(fn) =  commands.
% revised pab 21.11.2003
% - streamlined some code
% - updated help header
% - fixed bug in example  
% revised by es 25.05.00, error if frequencies are not equidistant or  do not
%                         start from zero 
% revised by es 23.05.00, call of freqtype, R.norm=S.norm  
% revised by pab 18.11.1999
%   - fixed a bug when S=S(f)
% revised by es 13.10.1999 
% revised by pab 23.08.1999

if ~isstruct(S)
  error('Incorrect input spectrum, see help datastructures')
end
  
if nargin<2||isempty(nr)
  nr=0; % number of derivatives
end

ftype = freqtype(S);
freq  = S.(ftype);
n     = length(freq);

if all(freq>0) % for .type='k2d', negative frequencies are allowed
  error('Spectrum does not start at zero frequency/wave number.\n Correct it with specinterp, for example.')
end
if any(abs(diff(freq,2))>1.0e-8)
  error('Not equidistant frequencies/wave numbers in spectrum.\n Correct it with specinterp, for example.')
end

if nargin<4||isempty(rate),
  rate=1; %interpolation rate
else
  if rate>16
    rate=16;
  end
  rate=2^nextpow2(rate);%make sure rate is a power of 2
end
if nargin<3||isempty(Nt),
  Nt = rate*(n-1);
else %check if Nt is ok
  Nt = min(Nt,rate*(n-1)); 
end

if numel(S.S)~= length(S.S) && (nargin>4 && Nx==0) && (nargin>5 && Ny==0)
  S = spec2spec(S,'freq');
  % If Nx=Ny=0 then a twodimensional spectrum gives no extra information
  % transform to type 'freq', and you do not have to call dspec2dcov
end
if numel(S.S)~=length(S.S)||(nargin>4 && Nx>0)||(nargin>5 && Ny>0)
  if nargin < 5
    Nx=[];Ny=[];dx=[];dy=[];
  end
  if nargin < 6
    Ny=[];
  end
  if nargin < 7
    dx=[];
  end
  if nargin < 8
    dy=[];
  end
  
  R = dspec2dcov(S,nr,Nt,rate,Nx,Ny,dx,dy);
  return % !!!!
end

if strcmpi(ftype,'k')
  vari='x';
else
  vari='t';
end

R      = createcov(nr,vari);
%R.R    = zeros(Nt+1,1);
R.tr   = S.tr;
R.h    = S.h;
R.norm = S.norm;
R.note = S.note;


%normalize spec so that sum(specn)/(n-1)=R(0)=var(X)
switch lower(ftype)
 case {'w','k'},
  w     = freq(:);
  dT    = pi/w(n);
  specn = S.S(:)*freq(n); %S.S(:)*pi/dT;
 case 'f',
  w     = 2*pi*freq(:);
  dT    = 1/(2*freq(n));  % sampling interval=1/Fs
  specn = S.S(:)*freq(n); %S.S(:)/(2*dT);
 otherwise
  error('unknown frequency type')
end

nfft = rate*2^nextpow2(2*n-2);

Rper = [specn; zeros(nfft-(2*n)+2,1) ; conj(specn(n-1:-1:2))]; % periodogram
t    = (0:Nt)'*dT*((2*n-2)/nfft);

%R   = setfield(R,vari,t);
R.(vari) = t;
r   = real(fft(Rper,nfft))/(2*n-2);
%r   = real(fft(Rper/(2*n-2),nfft));
R.R = r(1:Nt+1); 
if nr>0
  w         = [w ; zeros(nfft-2*n+2,1) ;-w(n-1:-1:2) ];
  fieldname = ['R' vari(ones(1,nr)) ];
  for ix=1:nr, 
    Rper = (-i*w.*Rper);
    r    = real(fft(Rper,nfft))/(2*n-2);
    R.(fieldname(1:ix+1)) = r(1:Nt+1);
    %R    = setfield(R,fieldname(1:ix+1),r(1:Nt+1));
    %eval(['R.R' vari(ones(1,ix)) '=r(1:(Nt+1));']);   
  end  
end






