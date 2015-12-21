function scaleflag=plotscaleflag(self,scale)
% DATA_1D/PLOTSCALEFLAG Return plotflag from plotscale
%
% CALL scaleflag = plotscaleflag(scale)
%
% scale    = string defining plotscale valid options are:
%           'linear', 'xlog', 'ylog', 'xylog', 'zlog', 'xzlog',
%           'yzlog',  'xyzlog' 
% scaleflag = integer defining plotscale.
%        0 'linear' : Linear scale on all axes.
%      100 'xlog'   : Log scale on x-axis.
%      200 'ylog'   : Log scale on y-axis.
%      300 'xylog'  : Log scale on xy-axis.
%      400 'zlog'   : Log scale on z-axis.
%      500 'xzlog'  : Log scale on xz-axis.
%      600 'yzlog'  : Log scale on yz-axis.
%      700 'xyzlog' : Log scale on xyz-axis.
%
%
% Example
% plotscaleflag(data_1d,'xyzlog')  % 700
% plotscaleflag(data_1d,'lin')     % 0
% plotscaleflag(data_1d,'ylog')    % 200
%
% See also data_1d/plotscale


% Note the ordering of validnames can not be changed without changing
% the order in functions dependent on this function

validScales = {'linear','xlog','ylog','xylog','zlog','xzlog','yzlog','xyzlog'};
if iscell(scale)
  scale = strvcat(scale{:});
end
  

if ischar(scale)
  logXscale = any(scale=='x',2);
  logYscale = any(scale=='y',2);
  logZscale = any(scale=='z',2);
%  scaleflag = 100*(logXscale + 10*(logYscale + 10*logZscale));
  scaleflag = 100*(logXscale +  2*(logYscale +  2*logZscale));
elseif isnumeric(scale)
  N = size(validScales,1);
  scaleflag = scale;
  scaleflag(scaleflag>N) = nan;
else
    error('Input must be character- or cell- arrays!')
end
  
%scaleflag = (scaleflag-1)*100;
return


