function Snew = specinterp(S,dt,Nmin,Nmax,method)
%SPECINTERP Interpolation and zero-padding of spectrum
%            to change Nyquist freq.
%
% CALL:  Snew = specinterp(S,dt,Nmin,Nmax,method)
%
%    Snew, S = spectrum structs (any type)
%         dt = wanted sampling interval (default as given by S, see spec2dt)
%              unit: [s] if frequency-spectrum, [m] if wave number spectrum  
%       Nmin = minimum number of frequencies. (default 0)
%       Nmax = minimum number of frequencies  (default 2^13+1 = 8193)
%      method = interpolation method (see interp1). (default 'pchip')
%              
% To be used before simulation (e.g. spec2sdat) or evaluation of covariance
% function (spec2cov) to directly get wanted sampling interval.
% The input spectrum is interpolated and padded with zeros to reach
% the right max-frequency, w(end)=pi/dt, f(end)=1/(2*dt), or k(end)=pi/dt.
% The objective is that output frequency grid should be at least as dense as
% the input grid, have equidistant spacing and length equal to 2^k+1 (>=Nmin).
% If the max frequency is changed, the number of points in the spectrum is maximized to
% 2^13+1. 
% NB! Also zero-padding down to zero freq, if S does not start there. 
%     If empty input dt, this is the only effect.
%  
% See also  spec2cov, spec2sdat, covinterp, spec2dt

% Tested on: Matlab 2013, 2015, 2016
% History:
% revided by GL May 2017 to new definition of fplot
% revised pab jan 2007
%  -added Nmin,Nmax, method, doInterpolate
%  -extra check on spacing using fplot
%  - changed default interpolation from linear to pchip.
% revised pab 12.10.2001
% -fixed a bug created 11.10.2001: ftype='k' now works OK
% revised pab 11.10.2001
% - added call to freqtype.m
% - fixed a bug: ftype=='f' now works correctly
% revised es 080600, revision of last revision, output matrix now OK
%  revised by es 22.05.00, output vectors columns, not rows
%  by es 13.01.2000, original ideas by sylvie, ir
  
Snew  = S;  
ftype = freqtype(S);
w     = S.(ftype);
n     = length(w);

%doInterpolate = 0;

if strcmp(ftype,'f') %ftype==f
  Cnf2dt = 0.5; % Nyquist to sampling interval factor
else % ftype == w og ftype == k
  Cnf2dt = pi;
end
wnOld = w(n);         % Old Nyquist frequency
dTold = Cnf2dt/wnOld; % sampling interval=1/Fs


if nargin<2||isempty(dt),  
  dt=dTold; 
end
if nargin<3||isempty(Nmin)
  Nmin = 0;
end
if nargin<4||isempty(Nmax)
  Nmax = 2^13+1;
end
if nargin<5||isempty(method)
  method = 'pchip';
end


% Find how many points that is needed
nfft   = 2^nextpow2(max(n-1,Nmin-1));
dttest = dTold*(n-1)/nfft;

while (dttest>dt) && (nfft<Nmax-1)
 nfft   = nfft*2;
 dttest = dTold*(n-1)/nfft;
end;
nfft = nfft+1;

wnNew         = Cnf2dt/dt; % New Nyquist frequency
dWn           = wnNew-wnOld;
doInterpolate = dWn>0 || w(1)>0 || (nfft~=n) || dt~=dTold || any(abs(diff(w,2))>1.0e-8); 

if doInterpolate>0
  if size(S.S,2)<2
    S1 = S.S.'; % if vector, make it a row to match matrix (np x nf)
  else
    S1 = S.S;
  end

  w  = w(:);
  %dw = w(n)-w(n-1);
  dw = min(diff(w));

  if dWn>0
    % add a zero just above old max-freq, and a zero at new max-freq
    % to get correct interpolation there
    Nz = 1  + (dWn>dw); % Number of zeros to add
    if Nz==2
      w = [w; wnOld+dw; wnNew];
    else
      w = [w; wnNew];
    end
    S1 = [S1 zeros(size(S1,1),Nz)];
  end
  if w(1)>0
    % add a zero at freq 0, and, if there is space, a zero just below min-freq
    Nz = 1 + (w(1)>dw); % Number of zeros to add
    if Nz == 2
      w=[0; w(1)-dw; w];
    else
      w=[0; w];
    end
    S1 = [zeros(size(S1,1),Nz) S1];
  end
  
  % Do a final check on spacing in order to check that the gridding is
  % sufficiently dense:
  np = size(S1,1);
  dwMin = realmax;
  %wnc = min(wnNew,wnOld-1e-5);
  wnc = wnNew;
  figure(99)
  for ix = 1:np
      if verLessThan('matlab','9.1')
          x = fplot(@(x)evalspec(x,ix),[0,wnc]);
      else
          HH = fplot(@(x)evalspec(x,ix),[0,wnc]);
          x = HH.XData;
          clear HH
      end
      dwMin = min(min(diff(x)),dwMin);
  end
  close(99)
  newNfft = 2^nextpow2(ceil(wnNew/dwMin))+1;
  if newNfft>nfft
    if (nfft<=2^15+1) && (newNfft>2^15+1)
      warning('WAFO:SPECINTERP','Spectrum matrix is very large (>33k). Memory problems may occur.')
    end
    nfft = newNfft;
  end
  
  Snew.(ftype) = linspace(0,wnNew,nfft).';
 
  Snew.S       = interp1(w,S1.',Snew.(ftype),method).';
  
end

if size(Snew.S,1)<2
  Snew.S=Snew.S.'; % if vector, make it a column
end

  function Sout = evalspec(win,thi)
    Sout = interp1(w,S1(thi,:),win,'pchip',0);
  end %
end % function specinterp


