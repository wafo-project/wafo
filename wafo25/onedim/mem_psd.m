function [ret_psd,ret_freq]=mem_psd(ar_coeffs,residual,freq,sample_f,method)
%
% [ret_psd,ret_freq] = mem_psd(ar_coeffs,residual,freq,sample_f,method)
%
% Calculates the power spectrum of the autoregressive filter
%
%                          M
%  x(n) = residual.e(n) + SUM ar_coeffs(k).x(n-k)
%                         k=1
%  where x(n) is the filter output and e(n) is white noise input.
%  If the "freq" argument is a vector (of frequencies) the spectrum is
%  calculated using the polynomial method and the "method" argument is
%  ignored.  For scalar "freq", an integer power of 2, or "method='FFT'",
%  causes the spectrum to be calculated by FFT.  Otherwise, the spectrum
%  is calculated as a polynomial.  Note that it is more computationally
%  efficient to use the FFT method if length of the filter is not much
%  smaller than the number of frequency values. The spectrum is scaled so
%  that spectral energy between zero frequency and the Nyquist frequency
%  is the same as the time-domain energy (i.e. mean square of the signal).
%
% Arguments:
%   ar_coeffs = [real vector] list of M=(order+1) autoregressive filter
%                N.B. the first element of "coeffs" is the zero-lag
%                coefficient, which always has a value coeffs(1)=1.
%
%   residual  = [real scalar] the moving-average coefficient of the AR
%                            filter.
%   freq      = [real vector] frequencies at which power spectral density
%                            is calculated
%              [integer scalar] number of frequency values (uniformly
%                      distributed from zero to the Nyquist frequency) at
%                      which spectral density is calculated. [default=256]
%   sample_f  = [real scalar] sampling frequency (Hertz) [default=1]
%   method    = [string] controls the method of evaluation if "freq"
%                       is a scalar ---
%              method="FFT":  use FFT to calculate power spectrum.
%              method="poly": use polynomial method
%              [default="poly" unless number of frequencies is an
%                         integer power of 2
%   freq, sample_f and method are optional.
%   freq and sample_f may be empty.
%
%  Returned values:
%     If return values are not required by the caller, the spectrum
%     is plotted and nothing is returned.
% ret_psd     = [real vector] power-spectrum estimate
% ret_freq    = [real vector] frequency values
%

% REFERENCES
% William H. Press and Saul A. Teukolsky and William T. Vetterling and
%               Brian P. Flannery",
% "Numerical recipes in C, The art of scientific computing", 2nd edition,
%    Cambridge University Press, 2002 --- Section 13.7.
% N.B. The algorithm in Press et al. expects prediction-filter coefficients.
%      mem_spect requires the whitening-filter coefficients.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2006 (C) Peter V. Lanspeary
%
% mem_spect.m is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the Free
% Software Foundation; either version 2, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
% (http://www.gnu.org/copyleft/gpl.html) for more details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% sanity checks
if ( nargin > 2 )
  user_freqs = isvector(freq) && length(freq)>1;
end
if ( nargin < 2 )
  error( 'mem_spect(ar_coeffs,residual,freq,sample_f,method) needs >=2 args');
elseif ( ~ isvector(ar_coeffs) || ~ isreal(ar_coeffs) || length(ar_coeffs)<2 )
  error( 'mem_spect: error: arg 1 (ar_coeffs) not real vector of length>=2.' );
elseif ( ~ isscalar(residual) || ~ isreal(residual) || residual <= 0 )
  error( 'mem_spect: error: arg 2 (residual) is not real scalar > 0' );
elseif ( nargin > 2 && ~isscalar(freq) && ~user_freqs && ~isempty(freq) )
  error( 'mem_spect: error: arg 3 (freq) is not vector or scalar.');
elseif ( nargin > 2 && ~user_freqs && ~isempty(freq) && ...
    ( ~isreal(freq) || fix(freq)~=freq || freq <= 2 || freq >= 1048576 ) )
  error( 'mem_spect: error: arg 3, (freq) is not integer >=2, <=1048576.)' );
elseif ( nargin > 2 && user_freqs && ~isempty(freq) && ...
    ( ~all(isreal(freq)) || any(freq<0) ) )
  error( 'mem_spect: error: arg 3, (freq) values must be real and >=0' );
elseif ( nargin > 3 && ~isempty(sample_f) && ...
    ( ~isscalar(sample_f) || ~isreal(sample_f) || sample_f<0 ) )
  error( 'mem_spect: error: arg 4, sample_f must be real scalar >0.' );
elseif ( nargin > 4 && ~ ischar(method) && ...
    ~ strcmp(method,'FFT') && ~ strcmp(method,'poly') )
  error( 'mem_spect: error: arg 5, (method) must be "FFT" or "poly".' );
else
  %% end of sanity checks
  %%
  %% define the frequencies
  if ( nargin < 3 || ( nargin >= 4 && isempty(freq) ) )
    freq = 256;
    user_freqs = 0;
  end
  if ( nargin < 4 || ( nargin >= 5 && isempty(sample_f) ) )
    sample_f = 1.0;
  end
  if ( user_freqs )
    len_freq = length(freq);
  else
    len_freq = freq;
    freq = (0:len_freq) * sample_f / 2 / len_freq;
  end
  %%
  %% decide which method to use
  is_power_of_2 = rem(log(len_freq),log(2))<10.*eps;
  force_FFT = nargin>4 && strcmp(method,'FFT');
  force_poly = nargin>4 && strcmp(method,'poly');
  use_FFT = ~user_freqs && ( ( ~ force_poly && is_power_of_2 ) || force_FFT);
  %%
  %% and do it
  len_coeffs = length(ar_coeffs);
  if ( use_FFT )
    fft_out = fft( [ ar_coeffs zeros(1,len_freq*2-len_coeffs ) ] );
    fft_out = fft_out(1:len_freq+1);
  else %% do not use FFT
    two_pi_i = (0+i) * 2 * pi / sample_f;
    fft_out = ar_coeffs * exp( two_pi_i * (1:len_coeffs)' * freq );
  end
  spectrum = 2.0 * residual / sample_f ./ ( fft_out .* conj(fft_out) );
  %%
  if ( nargout <= 0 )
    loglog(freq,spectrum);
  else
    ret_psd = spectrum;
    ret_freq = freq;
  end
end
end