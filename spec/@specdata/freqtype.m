function ftype=freqtype(S)
% FREQTYPE returns the frequency type of a Spectral density object.
%
% CALL: ftype = freqtype(S)
%
%  ftype = 'f' if frequency is given in Hz
%          'w' if frequency is given in rad/s
%          'k' if a wave number spectrum is given
%      S = spectral density struct
%
% See also  datastructures

% History:
%  pab 17.02.2007

ftype  = lower(S.freqtype);
