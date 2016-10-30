function ftype=freqtype(S)
% FREQTYPE returns the frequency type of a Spectral density struct.
%
% CALL: ftype = freqtype(S)
%
%  ftype = 'f' if frequency is given in Hz
%          'w' if frequency is given in rad/s
%          'k' if a wave number spectrum is given
%      S = spectral density struct
%
%  Example
%  S = demospec();
%  assert(freqtype(S), 'w')
%  S2 = ttspec(S);
%  assert(freqtype(S2), 'f')
%
% See also  datastructures

% History:
%  revised pab 17.02.2000
% - updated header
% revised pab 24.01.2000
%  - added check if field exist
% by pab 

names=fieldnames(S); 
ind=find(strcmp(names,'f')+strcmp(names,'w')+strcmp(names,'k')); %options are 'f' and 'w' and 'k' 
if isempty(ind)
  %ftype=[];
  error('WAFO:FREQTYPE','This is not a spectral density structure')
else
  ftype=lower(names{ind}); 
end