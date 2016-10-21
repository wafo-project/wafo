function [dt, wmdt,ftype] = spec2dt(S)
% SPEC2DT  Computes sampling interval from Nyquist frequency in spectrum
%
%  CALL : [dt,wmdt] = spec2dt(S)
%
%         dT = sampling interval, unit: 
%             [m] if wave number spectrum, 
%             [s] othrewise  
%       wmdt = 0.5 if wm have unit [Hz] 
%              pi  otherwise
%         S  = spectrum struct  
%
%  Let wm be maximum frequency/wave number in spectrum,
%  then dT=pi/wm if angular frequency, dT=1/(2*wm) if natural frequency (Hz)
%
% Example
%  S = jonswap;
%  dt = spec2dt(S)
% 
% See also

% Tested on Matlab 7.0
% revised pab 
% replaced getfield call with S.(ftype)
% Revised by es 25.05.00: help text + call to freqtype
  
ftype = freqtype(S);
if strcmp(ftype,'f') %ftype==f
  wmdt = 0.5; % Nyquist to sampling interval factor
else % ftype == w og ftype == k
  wmdt = pi;
end
wm = S.(ftype)(end); %Nyquist frequency
dt = wmdt/wm; % sampling interval=1/Fs


