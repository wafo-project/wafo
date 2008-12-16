function [coeffs,residual] = burg_filter(data,poles,stop_crit)
%BURG_FILTER Estimate 
%
% CALL [coeffs,residual] = burg_filter(data,poles,stop_crit)
%
% Calculates coefficients of a moving average (MA) whitening filter
% using the lattice-filter of Burg (1968) and also described by Kay
% and Marple (1981).  This filter reduces the data to white noise.
% Two other filters can be built from this whitening filter:
%   (1) a moving average linear prediction filter -- by removing the
%       first (zero-lag) coefficient and negating the others.
%   (2) an autoregressive generating (AR) filter -- by using "residual"
%       as the sole MA coefficient and using "coeffs" to supply the AR
%       coefficients of an ARMA filter.
% The power spectrum of the ARMA filter is an estimate of the maximum
% entropy power spectrum of the data.
%
% Arguments:
%   data      = [real vector] sampled data
%   poles     = [integer scalar] required number of poles of AR filter
%   stop_crit = [string] if 'FPE' or 'AIC', the FPE or AIC criterion
%                (Kay & Marple, 1981) are used to override the 'poles'
%                 argument so the filter does not grow too long.
%
% Returned values:
%   coeffs    = [real vector] list of M=(poles+1) moving-average filter
%                            coefficients; for data input x(n) and
%                            white noise output e(n), the filter is
%                        M
%                e(n) = SUM coeffs(k).x(n-k)
%                       k=0
%              N.B. the first element of "coeffs" is the zero-lag
%              coefficient, which always has a value coeffs(1)=1.
%
%   residual  = mean square of residual (white) noise from filter
%

% REFERENCES
% John Parker Burg (1968):
%   "A new analysis technique for time series data",
%   NATO advanced study Institute on Signal Processing with Emphasis on
%   Underwater Acoustics, Enschede, Netherlands, Aug. 12-23, 1968.
%
% Steven M. Kay and Stanley Lawrence Marple Jr. (1981):
%   "Spectrum analysis -- a modern perspective",
%   Proceedings of the IEEE, Vol 69, pp 1380-1419, Nov., 1981
%

% Please report bugs to pvl@...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2006 (C) Peter V. Lanspeary
%
% burg_filter.m is free software; you can redistribute it and/or modify it
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
%
% sanity checks
if ( nargin < 2 )
  error( 'burg_filter(data,poles): Need at 2 args.' );
elseif ( ~ isvector(data) || ~ isreal(data) || length(data) < 3 )
  error( 'burg_filter:error: arg 1 (data) must be real vector of length >3.' );
elseif ( ~isscalar(poles) || ~isreal(poles) || fix(poles)~=poles || poles<=0.5)
  error( 'burg_filter:error: arg 2 (poles) must be positive integer.' );
elseif ( floor(poles+0.5) > length(data)-2 )
  error( 'burg_filter:error: arg 2 (poles)+2 must be less than data length' );
elseif ( nargin>2 && ( ~ischar(stop_crit) || size(stop_crit,1)>1 ) )
  error( 'burg_filter:error: arg 3 (stop_crit) must be one string' );
else
%% end of sanity checks
%%
coeffs = zeros(1,poles);
%%
%% Storage of forward and backward prediction errors is a little tricky.
%% Because the forward error e(i) is always combined with the lagged
%% backward error b(i-1), e(1) and b(n) are never used, and therefore are
%% never stored.  Not storing unused data makes the calculation of the
%% reflection coefficient look much cleaner :)
data = reshape(data,1,[]);
N = size(data, 2);
forw_err = data(2:N);
back_err = data(1:N-1);
residual = sumsq(data)/N;
%%
%% Decide if FPE or AIC criteria will be applied to stop before
%% getting to the specified number of filter poles.
if ( nargin > 2 )
  use_FPE = strcmp(stop_crit,'FPE');
  use_AIC = strcmp(stop_crit,'AIC');
else
  use_FPE = 0;
  use_AIC = 0;
end
new_criterion = residual;
%old_criterion = 2 * new_criterion;
for k = 1:poles
  %%
  %%  reflection_coeff = -2 * E(e(i)*b(i-1)) / ( E(e(i)^2) + E(b(i-1)^2) )
  refl_coeff= -2 *forw_err * back_err' / (sumsq(forw_err) + sumsq(back_err));
  %%  Levinson-Durbin recursion for residual
  new_residual = residual * ( 1.0 - refl_coeff^2 );
  if ( k > 1 )
    %%
    %% Apply the FPE or AIC criterion and stop if the FPE or AIC is
    %% increasing rather than decreasing.
    %% Do it before we update the old filter "coeffs" and "residual".
    if ( use_FPE )
      old_criterion = new_criterion;
      new_criterion = new_residual * ( N + k + 1 ) / ( N - k - 1 );
      if ( new_criterion > old_criterion )
        break;
      end
    elseif ( use_AIC )
      old_criterion = new_criterion;
      new_criterion = log(new_residual) + 2 * ( k + 1 ) / N;
      if ( new_criterion > old_criterion )
        break;
      end
    end
    %% Update filter "coeffs" and "residual".
    %% Use Levinson-Durbin recursion formula.
    coeffs = [ prev_coeffs+refl_coeff.*prev_coeffs(k-1:-1:1), refl_coeff ];
    residual = new_residual;
  else
    coeffs = refl_coeff;
    residual = new_residual;
  end
  if ( k < poles )
    prev_coeffs = coeffs;
    %%  calculate new prediction errors (by recursion):
    %%  e(p,i) = e(p-1,i)   + k * b(p-1,i-1)  i=2,3,...n
    %%  b(p,i) = b(p-1,i-1) + k * e(p-1,i)    i=2,3,...n
    %%  remember e(p,1) is not stored, so don't calculate it; make e(p,2)
    %%  the first element in forw_err.  b(p,n) isn't calculated either.
    nn = length(forw_err);
    forw_new = forw_err(2:nn)   + refl_coeff .* back_err(2:nn);
    back_err = back_err(1:nn-1) + refl_coeff .* forw_err(1:nn-1);
    forw_err = forw_new;
  end
end
%% end of for loop
coeffs = [1 coeffs];
end
end