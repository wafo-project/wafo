function scale=plotscale(self,plotflag)
% DATA_1D/PLOTSCALE Return plotscale from plotflag
%
% CALL scale = plotscale(plotflag)
%
% plotflag = integer defining plotscale.
%   Let scaleId = floor(plotflag/100). 
%   If scaleId < 8 then:
%      0 'linear' : Linear scale on all axes.
%      1 'xlog'   : Log scale on x-axis.
%      2 'ylog'   : Log scale on y-axis.
%      3 'xylog'  : Log scale on xy-axis.
%      4 'zlog'   : Log scale on z-axis.
%      5 'xzlog'  : Log scale on xz-axis.
%      6 'yzlog'  : Log scale on yz-axis.
%      7 'xyzlog' : Log scale on xyz-axis.
%  otherwise
%   if (mod(scaleId,10)>0)            : Log scale on x-axis.
%   if (mod(floor(scaleId/10),10)>0)  : Log scale on y-axis.
%   if (mod(floor(scaleId/100),10)>0) : Log scale on z-axis.
%
% scale    = string defining plotscale valid options are:
%       'linear', 'xlog', 'ylog', 'xylog', 'zlog', 'xzlog',
%       'yzlog',  'xyzlog' 
%
% Example
% plotscale(data_1d,100)  % xlog
% plotscale(data_1d,200)  % ylog
% plotscale(data_1d,1000) % ylog
%
% See also data_1d/plotscaleflag
scaleId = floor(plotflag/100);
if scaleId<8
  scaleId = scaleId+1;
else
  logXscaleId = (mod(scaleId,10)>0);

  logYscaleId = (mod(floor(scaleId/10),10)>0)*2;
  logZscaleId = (mod(floor(scaleId/100),10)>0)*4;
  scaleId = logYscaleId +logXscaleId+logZscaleId +1;
end
scales = {'linear','xlog','ylog','xylog','zlog','xzlog','yzlog','xyzlog'};

scale = scales{scaleId};