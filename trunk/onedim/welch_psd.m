function [varargout] = welch_psd(varargin)
%WELCH_PSD  Estimate power spectral density using Welch's method
%
% CALL: [Pxx, f, options, Pci] = welch_psd(x,options,...)
%
% Pxx,Pci = Power spectral density and confidence interval, respectively.
% f       = frequency vector 
% x       =  vector/matrix of samples
% options = struct defining spectral estimation. Valid parameters are:
%       nfft    : size of fourier transform window (default min(256,length(x)))
%       Fs      : sampling frequency,  (default 2 Hz)
%       window  : vector, string or function handle defining the window function
%                 Note: window length can be specified instead, in which case
%                 window = hanning(length)    (default 'hanning')
%       overlap : overlap ratio or number of overlap samples between sections,
%                 0    = means no overlap, 
%                 0.5  = means 50% overlap (default)
%                 0.75 = means 75% overlap 
%                
%       p       : confidence level for confidence interval (default=0)
%                 P must be between 0 and 1; no confidence interval will be computed for P == 0.
%       range   : 'twosided' : return all frequencies           (default if x is complex) 
%                 'onesided' : return half of the frequencies   (default if x is real)
%       dflag   : 'none','mean', 'linear', 'ma' (moving average) (default='none')
%                 remove trends from the data sections before computing spectral estimates
%
% WELCH_PSD estimate power spectrum of a stationary signal. This chops the 
% signal into overlapping sections, windows each section and applies a Fourier
% transform to determine its frequency components. The magnitudes of these sections
% are then averaged to produce the estimate Pxx.
% The confidence interval around the estimate is returned in Pci.
%
%
% Examples
%    [b,a] = cheby1(4,3,[0.2, 0.4]);    % define noise colour
%    y = filter(b,a,randn(2^12,1));
%    welch_psd(y);                      % estimate noise colour
%    opt = welch_psd('defaults');
%    opt.nfft = 128;
%    opt.Fs = 3;
%    welch_psd(y,opt);
%    welch_psd(y,'nfft',128,'Fs',3)    % alternatively
%
% See also welch_cpsd

% History
% revised pab 2007
% -Made it conform to wafo style.
% -changed name from pwelch to welch_psd
% -Made it possible to compute PSD for matrices as well.
% -changed input options
% 2001-04-02 Paul Kienzle
%    * return nfft/2+1 elements rather than nfft/2 for even nfft.
%    * use more accurate (and faster) computation of confidence intervals

% TODO % does not work correctly yet for matrix inputs
% TODO: Should be extended to accept a vector of frequencies at which to
% TODO:    evaluate the fourier transform (via filterbank or chirp
% TODO:    z-transform).
% TODO: What should happen with the final window when it isn't full?
% TODO:    currently I dump it, but I should probably zero pad and add
% TODO:    it in.
% TODO: Consider returning the whole gamit of Pxx, Pyy, Pxy, Cxy, Txy
% TODO:    as well as confidence intervals for each;  if users tend
% TODO:    only to use one of these don't bother, but if they ever need
% TODO:    more than one, then it's free.  Alternatively, break out the
% TODO:    compute engine into a function that the user can call directly.
% TODO: Check if Cxy, Txy are computed frame-by-frame or on the average
% TODO:    of the frames.  SpcTools and I do it on the average, 
% TODO:    wdkirby@ix.netcom.com (1998-04-29 octave-sources) computes 
% TODO:    them frame-by-frame.



% Copyright (C) 1999-2001 Paul Kienzle
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


% sort out parameters
if nargin < 1,
  disp('[Pxx, f] = welch_psd(x,options)');
  return
end

options = struct('nfft',[],'window','hanning','overlap',0.5,'Fs',2,'range',[],'dflag','none','p',0);

if nargin==1 && strcmpi(varargin{1},'defaults')
  varargout{1} = options;
  return
end


if isa(varargin{1},'cell')
  % Secret call
  %  [P,f] = welch_psd({x,y},method,options....)
  try
    x = varargin{1}{1};
    y = varargin{1}{2};
    method = varargin{2}; % Determine if we are called as cpsd, cohere or tfe
    switch method
      case 'cpsd'
        ftype = 2;
      case 'cohere'
        ftype = 3;
      case 'tfe'
        ftype = 4;
      otherwise
        error('WAFO:SPEC:WELCH','Unknown option: %s',method)
    end
  catch
    error('Input must be given as [] = welch_psd({x,y},method,options...)')
  end
  first = 3;
else
  %method = 'psd';
  ftype = 1;
  first = 2;
  x=varargin{1};
  y = [];
end

[Nx,Mx] = size(x);
if Nx == 1,
  x = x.';
  [Nx,Mx] = size(x);
end

if ~isempty(y)
  [Ny,My] = size(y);
  if Ny == 1,
    y = y.';
    [Ny,My] = size(x);
  end
end
if (Nx <= 1 || Mx == 0) || ...
    (~isempty(y) && (Ny <= 1 && My < 1))
  error ([method, ' data must be a vector or matrix']);
end

if ~isempty(y) && (Nx~=Ny || Mx~=My)
  error ([method, ' x and y vectors/matrices must be of same size!']);
end

if nargin>=first
  options = parseoptions(options,varargin{first:nargin});
end

% Fill in defaults for arguments that aren't specified
if isempty(options.nfft), 
  options.nfft = min(256, length(x)); 
end

switch class(options.window)
  case {'function_handle', 'char','inline'}
    window = feval(options.window,options.nfft);
  otherwise % numeric input
    if length(options.window) ==1
      % if only the window length is given, generate hanning window
      window = hanning(options.window);
    else
      window = options.window;
    end
end

if size(window,1)==1, 
  window = window.';
end

if options.overlap<1
  noverlap = fix(length(window)*options.overlap);
else
  noverlap = fix(options.overlap);
end
window = window(:,ones(1,Mx));
  
if isempty(options.range),
  if ~isreal(x) || (~isempty(y) && ~isreal(y))
    options.range = 'twosided';
  else
    options.range = 'onesided';
  end
end

if strncmpi(options.dflag,'none',1)
  options.dflag = 'none';
  trend=-1;
elseif strncmpi( options.dflag,'constant',1) || strncmpi( options.dflag,'mean',2)
  options.dflag = 'mean';
  trend = 0;
%   x = detrend(x,trend);
%   if ~isempty(y)
%     y = detrend(y,trend);
%   end

elseif strncmpi( options.dflag,'ma',2)
  trend = 0;
  options.dflag = 'ma';
%   x = detrendma(x,2*options.nfft);% remove moving average
%   if ~isempty(y)
%      y = detrendma(y,2*options.nfft);
%    end

  
else %if strncmpi( options.dflag,'linear',1)
  options.dflag = 'linear';
  trend = 1;
%    x = detrend(x,trend);
%    if ~isempty(y)
%      y = detrend(y,trend);
%    end
end
 

  if (ftype > 2 && options.p > 0)
    error([ method, ' can''t compute confidence intervals' ]);
  elseif (options.p < 0 || options.p > 1)
    error([ method, ' confidence interval must be between 0 and 1' ]);
  end

  % Normalize the window
  window = window / norm(window);

  % compute window offsets
  win_size = length(window);
  if (win_size > options.nfft)
    options.nfft = win_size;
    warning('WAFO:SPEC:WELCH','%s fft size adjusted to %d', method, win_size);
  end
  nfft = options.nfft;
  if win_size<=noverlap
    error('WAFO:SPEC:WELCH','%s overlap size must be smaller than window length %d', method, win_size);
  end
  step = win_size - noverlap;

  % Determine which correlations to compute
  Pxx = [];
  Pyy = [];
  Pxy = [];
  if ftype~=2, Pxx = zeros(nfft,Mx); end % Not needed for csd
  if ftype==3, Pyy = zeros(nfft,Mx); end % Only needed for cohere
  if ftype~=1, Pxy = zeros(nfft,Mx); end % Not needed for psd

  % Average the slices
  offset = 1:step:length(x)-win_size+1;
  N = length(offset);
  for i=1:N
    a = x(offset(i):(offset(i)+win_size-1),:);
    if trend>=0, 
      a=detrend(a,trend); 
    end
    a = fft(postpad(a.*window, nfft));
    if ~isempty(Pxx), 
      Pxx = Pxx + a.*conj(a);
    end
    if ~isempty(Pxy)
      b = y(offset(i):offset(i)+win_size-1,:);
      if trend>=0,
        b=detrend(b,trend); 
      end
      b = fft(postpad(b.*window, nfft));
      Pxy = Pxy + a .*conj(b);
      if ~isempty(Pyy),
        Pyy = Pyy + b.*conj(b); 
      end
    end
  end
  if (ftype <= 2)
    % the factors of N cancel when computing cohere and tfe
    if ~isempty(Pxx), Pxx = Pxx / N; end
    if ~isempty(Pxy), Pxy = Pxy / N; end
    if ~isempty(Pyy), Pyy = Pyy / N; end
  end

  % Compute confidence intervals
  if options.p > 0, 
    Pci = zeros(nfft,Mx); 
  end
  if (options.p > 0 && N > 1)
    if ftype>2
      error([method, ': internal error -- shouldn''t compute Pci']); 
    end

    % c.i. = mean +/- dev
    % dev = z_ci*std/sqrt(n)
    % std = sqrt(sumsq(P-mean(P))/(N-1))
    % z_ci = normal_inv( 1-(1-ci)/2 ) = normal_inv( (1+ci)/2 );
    % normal_inv(x) = sqrt(2) * erfinv(2*x-1)
    %    => z_ci = sqrt(2)*erfinv(2*(1+ci)/2-1) = sqrt(2)*erfinv(ci)
    for i=1:N
      a=x(offset(i):offset(i)+win_size-1,:);
      if trend>=0, 
        a=detrend(a,trend);
      end
      a=fft(postpad(a.*window, nfft));
      if ftype == 1 % psd
      	P = a.*conj(a) - Pxx;
      	Pci = Pci + P.*conj(P);
      else          % csd
      	b = y(offset(i):offset(i)+win_size-1,:);
      	if trend>=0,
          b=detrend(b,trend);
        end
      	b = fft(postpad(b.*window, nfft));
      	P = a.*conj(b) - Pxy;
      	Pci = Pci + P.*conj(P);
      end
    end
      
    Pci = ( erfinv(options.p) * sqrt( 2/N/(N-1) ) ) * sqrt ( Pci );
  end



Fs = options.Fs;
f = (0:nfft).'*Fs/nfft;

  switch (ftype)
    case 1, % psd      
     [P,f] = onesidedOrTwosided(Pxx/Fs,f,options); 
    case 2, % cpsd
      [P,f] = onesidedOrTwosided(Pxy/Fs,f,options);
    case 3, % cohere
      [P,f] = onesidedOrTwosided(Pxy.*conj(Pxy)./Pxx./Pyy,f,options);
    case 4, % tfe
      [P,f] = onesidedOrTwosided(Pxy./Pxx,f,options);
  end

  if options.p > 0,
    [Pci,f] = onesidedOrTwosided(Pci/Fs,f,options);
    Pci = [ P - Pci, P + Pci ];
  end
  
  
  % Plot if there is no 
  if nargout==0,
    %   unwind_protect
    if options.p>0
      plot(f, 10*log10(abs(P)), 'b-',f, 10*log10(abs(Pci)), 'r-');
    else
      plot(f, 10*log10(abs(P)), 'b-');
    end
    if Fs==2
      xlabel('Frequency (rad/pi)');
    else
      xlabel('Frequency (Hz)');
    end
    if ftype==1
      title ('Welch''s Spectral Estimate Pxx/Fs');
      ytext='Power Spectral Density';
    elseif ftype==2
      title ('Cross Spectral Estimate Pxy');
      ytext='Cross Spectral Density';
    elseif ftype==3
      title ('Coherence Function Estimate |Pxy|^2/(PxxPyy)');
      ytext='Coherence ';
    else
      title ('Transfer Function Estimate Pxy/Pxx');
      ytext='Transfer';
    end
    %if use_dB,
      ylabel(strcat(ytext, ' (dB)'));
%     else
%       ylabel(ytext);
%     end
    grid('on');
    
%     unwind_protect_cleanup
%       grid('off');
%       title('');
%       xlabel('');
%       ylabel('');
%     end_unwind_protect
  end
	   
  if nargout>=1, 
    varargout{1} = P; 
    if nargout>1
      varargout{2} = f;
      if nargout>2
        varargout{3} = options;
        if nargout>3 && options.p>0,
          varargout{4} = Pci;
        end
      end
    end   
  end
end

function [P,f] = onesidedOrTwosided(P,f,options)
nfft = options.nfft;
if strncmpi(options.range,'onesided',1)
  if rem(nfft,2)==1 % Odd
    n = (nfft+1)/2;
    P = [P(1,:); P(2:n,:)*2];
  else   % Even
    n = nfft/2 + 1; % include DC and Nyquist
    P = [P(1,:); P(2:n-1,:)*2; P(n,:)];
  end
  f = f(1:n);
else % twosided
   %donothing
   % n = nfft;
end

end % onesidedOrTwosided




%~demo
% Fs=8000;
% [b,a] = cheby1(4,3,2*[500, 1000]/Fs);    % define spectral envelope
% s=0.05*randn(2^11,1);                    % define noise
% idx=fix(1:Fs/70:length(s))'; 
% s(idx)=s(idx)+ones(size(idx));           % add 70 Hz excitation
% x=filter(b,a,s);                         % generate signal
%
% figure(1); subplot(221); 
% text(0,0.9,'basic estimate','Units','Normalized'); 
% pwelch(x',[],Fs); text;   % slip in a test for row vs. column vector
% subplot(222); 
% text(0,0.9,'nfft=1024 instead of 256','Units','Normalized'); 
% pwelch(x,1024); text;
% subplot(223); 
% text(0,0.9,'boxcar instead of hanning','Units','Normalized');
% pwelch(x,[],[],boxcar(256)); text;
% subplot(224); 
% text(0,0.9,'no overlap','Units','Normalized'); 
% pwelch(x,[],[],[],0); text;
%
% figure(2); subplot(121);
% text(0,0.9,'magnitude units, whole range','Units','Normalized'); 
% pwelch(x,'whole','squared'); text;
% subplot(122);
% text(0,0.9,'90% confidence intervals','Units','Normalized'); 
% pwelch(x,[],[],[],[],0.9); text;
% oneplot();
% %----------------------------------------------------------
% % plots should show a chebyshev bandpass filter shape
