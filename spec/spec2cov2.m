function R = spec2cov2(S,nr,Nt,dt)
%SPEC2COV2 Computes covariance function and its derivatives, alternative version
% 
%    CALL:  R = spec2cov2(S,nr,Nt,dt);
% 
%  R    = [R0, R1,...Rnr] matrix with autocovariance and its
%          derivatives, i.e., Ri (i=1:nr) are column vectors with
%          the 1'st to nr'th derivatives of R0.  size Nt+1 x Nr+1
%  S    = a spectral density structure (See datastructures)
%  Nr   = number of derivatives in output, Nr<=4          (default 0)
%  Nt   = number in time grid, i.e., number of time-lags. 
%         (default rate*(n-1)) where rate = round(1/(2*f(end)*dt)) or 
%          rate = round(pi/(w(n)*dt)) depending on S.
%  dt   = time spacing for R
%  
%   NB! This routine requires that the spectrum grid is equidistant
%       starting from zero frequency.
% Example:    
%      S = jonswap;
%      dt = 0.1;  
%      R = spec2cov2(S,3,256,dt);
%
% See also  spec2cov, specinterp, datastructures  
  
% History
% revised pab Sept 2005
%  - rate is sometimes zero -> spec2cov2 crash: made an ad-hoc solution.
% by pab 24.11.2003
% based on code from spec2XXpdf programmes
  
error(nargchk(1,4,nargin))
if ~isstruct(S)
  error('Incorrect input spectrum, see help datastructures')
end
  
if nargin<2||isempty(nr)
  nr=0; % number of derivatives
end
ftype = freqtype(S); %options are 'f' and 'w' and 'k'
freq  = S.(ftype);
n     = length(freq);

if nargin<4||isempty(dt),
  if strcmpi(ftype,'f')
    dt=1/(2*freq(n)); % sampling interval=1/Fs
  else
    dt=pi/(freq(n));
  end
end
if strcmpi(ftype,'f')
  rate = round(1/(2*freq(n)*dt)); % sampling interval=1/Fs
else
  rate = round(pi/(freq(n)*dt));
end
% TODO % rate sometimes is zero -> spec2cov2 sometimes crash
rate = max(rate,1);

if nargin<3||isempty(Nt),
  Nt = rate*(n-1);
else %check if Nt is ok
  Nt = min(Nt,rate*(n-1));
end

checkdt = 1.2*min(diff(freq))/2/pi;
switch ftype
 case {'w'},
  vari = 't'; 
 case {'f'}, 
  vari    = 't'; 
  checkdt = checkdt*2*pi;
 case {'k'},
  vari = 'x';
end
msg1 = sprintf(['The step dt = %g in computation of the density is too' ...
		' small.'],dt);
msg2 = sprintf(['The step dt = %g step is small, may cause numerical' ...
		' inaccuracies.'],dt);

if (checkdt < 2^-16/dt)
  disp(msg1)
  disp('The computed covariance (by FFT(2^K)) may differ from the theoretical.')
  disp('Solution:')
  error('use larger dt or sparser grid for spectrum.')
end

% Calculating covariances
%~~~~~~~~~~~~~~~~~~~~~~~~
S2 = specinterp(S,dt);
R2 = spec2cov(S2,nr,Nt,1,0,0);
R  = zeros(Nt+1,nr+1);
fieldname = ['R'  vari(ones(1,nr))];
idx = 1:Nt+1;
for ix = 1:nr+1
  fn = fieldname(1:ix);
  R(:,ix) = R2.(fn)(idx);
end

EPS0 = 0.0001;
cc  = R(1,1)-R(2,1)*(R(2,1)/R(1,1));
if Nt+1>=5
  %cc1=R(1,1)-R(3,1)*(R(3,1)/R(1,1))+R(3,2)*(R(3,2)/R(1,3));
  %cc3=R(1,1)-R(5,1)*(R(5,1)/R(1,1))+R(5,2)*(R(5,2)/R(1,3));
  
  cc2 = R(1,1)-R(5,1)*(R(5,1)/R(1,1));
  if (cc2<EPS0)
    warning('WAFO:SPEC2COV2',msg1)
  end
end
if (cc<EPS0)
  disp(msg2)
end
return