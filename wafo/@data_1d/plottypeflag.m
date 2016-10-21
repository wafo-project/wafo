function flag=plottypeflag(self,plottype)
% DATA_1D/PLOTTYEPEFLAG Return plotflag from plottype
%
% CALL plotflag = plottypeflag(plottype)
%
% plottype = string defining plottype valid options are:
%           'plot', 'stairs', 'stem', 'errorbar', 'bar', 'area'
% plotflag = integer defining plottype.
%      1 'plot'    : Linear plot.
%      2 'stairs'  : Stairs plot
%      3 'stem'    : Stem plot
%      4 'errorbar': Errorbar plot
%      5 'bar'     : bar graph
%      6 'area'    : area plot
%
%
% Examples
% plottypeflag(data_1d,'stairs')  % 2
% plottypeflag(data_1d,'plot')    % 1
% plottypeflag(data_1d,'bar')     % 5
% plottypeflag(data_1d,'ba')      % nan
%
% See also data_1d/plottype


% Note the ordering of validnames can not be changed without changing
% the order in functions dependent on this function

validNames = {'plot','stairs','stem','errorbar','bar','area'};
if iscell(plottype)
  plottype = strvcat(plottype{:});
end
  

if ischar(plottype)
   N1 = size(plottype,1);
    flag = repmat(nan,N1,1);
    
    for ix = 1:N1
      tmp =    strmatch(lower(deblank(plottype(ix,:))),validNames,'exact');
      if ~isempty(tmp)
        flag(ix) = tmp;
      end
    end
elseif isnumeric(scale)
  N = size(validNames,1);
  flag = plottype;
  flag(flag>N) = nan;
else
    error('Input must be character- or cell- arrays!')
end
  

return


