function [psd,freq]=welch_psd2(signal,window,Fs,overlap,padding,qual)
%WELCH_PSD2 Estimate power spectrum using Welch's method alternative 
%
% CALL [psd,freq]=welch_psd2(signal,window,Fs,overlap,padding,qual)
%
% signal   = [real vector] time-series signal
% window   = [real vector] of window-function values between 0 and 1; the
%                  signal segment has the same length as the window, or
%           [integer scalar] length of each signal segment (and
%                 unpadded FFT);
%           default value is sqrt(length(x)) rounded down to the
%           nearest integer power of 2, and forces qual="sloppy"
% Fs      = [real scalar] sampling frequency (Hertz); default=1.0
% overlap = [real scalar] segment overlap factor
%            0 <= overlap < 1, default is zero, 0.5 is "industry standard"
% padding = [integer scalar] number of samples of zero padding (per FFT)
%                 default is zero
% qual    = [string] Quality specifier
%              "sloppy"  FFT length is rounded up to the nearest integer
%                        power of 2, segment length is rounded up unless
%                        the "window" arg is specified as a vector
%                        FFT length is adjusted after addition of padding
%          default is to use exactly the segment lengths and padding
%          lengths (hence FFT length) specified in argument list
%
% WELCH_PSD estimate power spectrum of a time-series signal by the periodogram
% (FFT) method.  The "signal" is divided into segments of length
% equal to the length of "window", and each segment is multiplied by
% "window" before (optional) zero-padding and calculating its FFT.
% The power spectrum is the mean square of the FFTs, scaled so that area
% under the power spectrum is the same as the mean square of the signal.
%
% All but the first argument are optional.
% Any but the first and last may be empty.
%
% Returned values:
%     If return values are not required by the caller, the spectrum
%     is plotted and nothing is returned.
% psd  = [real vector] power-spectrum estimate
% freq = [real vector] frequency values
%         the length of both returned vectors is [length_of_FFT]/2+1
% 
% Example
% x = load('sea.dat');
% Fs = 1/diff(x(1:2,1));
% [Si,fi] = welch_psd2(x(:,2),[],Fs);
% plot(fi, Si);
%
% close all;
%
% See also dat2spec

%
% REFERENCE
%  Peter D. Welch (June 1967):
%   The use of fast Fourier transform for the estimation of power spectra:
%   a method based on time averaging over short, modified periodograms.
%   IEEE Transactions on Audio Electroacoustics, VOl AU-15(6), pp 70-73
%

% History
% revised pab april 2007
% -renamed from welch_psd to welch_psd2
% -updated help header to conform to the wafo-style


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2006 (C) Peter V. Lanspeary
%
% welch_psd.m is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Soft-
% ware Foundation; either version 2, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
% (http://www.gnu.org/copyleft/gpl.html) for more details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

if ( nargin >0 )
  x_len = length(signal);
  is_win = 0;
  if ( nargin > 1 )
    if ( isscalar(window) )
      is_win = 1;
    elseif ( isvector(window) )
      is_win = length(window);
    end
  end
end

%%
%% SANITY CHECKS
if ( nargin <= 0 )
  error( 'welch_psd(signal,window,overlap,padding): Need at least 1 arg.' );
elseif ( ~isvector(signal) || ~ isreal(signal) )
  error( 'welch_psd: error: arg 1 (signal) must be real vector.' );
elseif ( nargin > 1 && ~isempty(window) && ~is_win )
  error( 'welch_psd: error: arg 2 must be window vector or segment length.' );
elseif ( nargin > 1 && is_win==1 && ( ~isreal(window) || ...
    fix(window)~=window || window<4 || x_len<window ) )
  error( 'welch_psd: error: arg 2, window not integer, >4 & <=length(data).' );
elseif ( nargin > 1 && is_win>1 && ( ~all(isreal(window)) || any(window<0) ) )
  error( 'welch_psd: error: arg 2, window values must be real and >=0' );
elseif ( nargin > 2 && ~isempty(Fs) && ...
    ( ~isscalar(Fs) || ~isreal(Fs) || Fs<0 ) )
  error( 'welch_psd: error: arg 3, Fs must be real scalar >0.' );
elseif ( nargin > 3 && ~isempty(overlap) && ...
    ( ~isscalar(overlap) || ~isreal(overlap) || overlap<0 || overlap>=1 ) )
  error( 'welch_psd: error: arg 4, overlap not between 0 and 1' );
elseif ( nargin > 4 && ~isempty(padding) && ...
    ( ~isscalar(padding) || ~isreal(padding) || ...
    fix(padding)~=padding || padding<0 ) )
  error( 'welch_psd: error: arg 5, padding not integer >=0' );
elseif ( nargin > 5 && ( ~ischar(qual) || size(qual,1)>1 ) )
  error( 'welch_psd: error: arg 6, qual must be one string' );
else
  %% end of sanity checks
  %%
  log_two = log(2);
  nearly_one = 0.99999999999;
  is_sloppy = ~is_win || ( is_win==1 && nargin>5 && strcmp(qual,'sloppy') );
  %%
  %% calculate/adjust segment length
  if ( ~is_win )
    seg_len = 2 ^ ceil( log( sqrt(x_len) ) * nearly_one / log_two );
    window = 1;
    win_meansq = 1;
  elseif ( is_win==1 )
    seg_len = window;
    window = 1;
    win_meansq = 1;
  else
    window = reshape(window,1,[]);
    seg_len = length(window);
    win_meansq = sumsq(window) / seg_len;
  end
  %%
  if ( nargin<3 || ( nargin>=3 && isempty(Fs) ) )
    Fs = 1;
  end
  if ( nargin<4 || ( nargin>=4 && isempty(overlap) ) )
    overlap = 0;
  end
  if ( nargin<5 || ( nargin>=5 && isempty(padding) ) )
    padding = 0;
  end
  if ( x_len < seg_len )
    error( 'welch_psd: error: Signal is shorter than segment/window' );
  end
  %%
  %% calculate FFT length
  fft_len = seg_len + padding;
  if ( is_sloppy )
    fft_len = 2 ^ ceil( log( fft_len ) * nearly_one / log_two );
  end
  %%
  %% AVERAGE THE SIGNAL
  %%
  %% Use the same signal data as the periodogram to take into account the
  %% overlap and any unused data.
  %% N.B. that removing the mean for ALL the data makes the spectral energy
  %% a bit too high, removing a "local/segment" mean value from each segment
  %% makes spectral energy too low. What is the CORRECT method? IDNK.
  %%
  overlap = fix(seg_len * overlap);
  avg_signal = 0;
  n_fft = 0;
  for start_seg = (1:seg_len-overlap:x_len-seg_len+1)
    avg_signal = avg_signal + sum(signal(start_seg:start_seg+seg_len-1));
    n_fft = n_fft +1;
  end
  avg_signal = avg_signal / n_fft / seg_len;
  %%
  %% CALCULATE PERIODOGRAM
  %%
  %% When computing the FFT of real data, it is possible to increase the
  %% speed of the algorithm by
  %%   x = signal(start_seg:2:end_seg) + i * signal(start_seg+1:2:end_seg);
  %% and then unscrambling the vector returned by fft().  The length of
  %% the decreaes by a factor of 2.  However, it is quite possible that
  %% fft()  recognises real input and uses this (or a similar) trick anyway.
  %%
  psd_len = floor(fft_len/2)+1;
  psd = zeros(1,psd_len);
  x = zeros(1,fft_len);
  n_fft = 0;
  for start_seg = (1:seg_len-overlap:x_len-seg_len+1)
    x(1:seg_len) = window .* ...
    (signal(start_seg:start_seg+seg_len-1) - avg_signal );
    X = fft(x);
    n_fft = n_fft +1;
    psd   = psd + abs( X(1:psd_len) ).^2;
  end
  psd = psd / ( n_fft * seg_len * Fs * win_meansq / 2 );
  freq = (0:psd_len-1) * ( Fs / fft_len );
  if ( nargout <= 0 )
    loglog(freq,psd);
  %else
   % ret_psd = psd;
   % ret_freq = freq;
  end
end
end
