function type1=plottype(self,flag)
% DATA_1D/PLOTTYPE Return plottype from plotflag
%
% CALL plottype = plottype(plotflag)
%
% plottype = string defining plottype valid options are:
%           'linear', 'stairs', 'stem', 'errorbar', 'bar', 'area'
% plotflag = integer defining plottype.
%      1 'plot'  : Linear plot.
%      2 'stairs'  : Stairs plot
%      3 'stem'    : Stem plot
%      4 'errorbar': Errorbar plot
%      5 'bar'     : bar graph
%      6 'area'    : area plot
%
%
% Examples
% plottype(data_1d, 2)  % 'stairs'
% plottype(data_1d,2:3)  % {'stairs','stem'}
%
% See also data_1d/plottypeflag


% Note the ordering of validnames can not be changed without changing
% the order in functions dependent on this function

validNames = {'plot','stairs','stem','errorbar','bar','area',nan};
if isnumeric(flag)
  N = size(validNames,1);
  
  flag(flag<N | N<flag) = N+1;
  if length(flag)>1
    type1 = validNames(flag);
  else
    type1 = validNames{flag};
  end
else
    error('Input must be numeric!')
end
  

return


