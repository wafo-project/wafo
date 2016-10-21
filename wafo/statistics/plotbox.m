%PLOTBOX   Plot box-and-whisker diagram
%
% CALL : s = plotbox(data,options);
%
%       s = matrix with one column for each dataset as follows:
%           1  minimum
%           2  1st quartile
%           3  2nd quartile (median)
%           4  3rd quartile
%           5  maximum
%           6  lower confidence limit for median
%           7  upper confidence limit for median
%    data = matrix with one column for each dataset, or data is a cell
%           vector with one cell for each dataset.
% optoins = name-value pairs or structure with fieldnames 
%     notched = 0  rectangular box plot (default)
%               a  0<a<1  notch of the specified depth
%               1  notched box plot
%      symbol = symbol(s) for outlier values, i.e., points outside the 
%               whisker range. (default '+o')
%               symbol(1) is used for points between 1 and 2 times MWR
%               symbol(2) is used for points larger than 2 times MWR
%    vertical = 0, makes the boxes horizontal
%               1, makes the boxes vertical (default)
%    width    = scalar or text defining width of the boxes. Options are
%               'fixed', 'auto'. (default 'auto')     
%  maxwhisker = defines the length of the whiskers as a function of the
%               interquartile range (IQR) (default = 1.5) 
%               Note: MWR = maxwhisker*IQR.
%
% The box plot is a graphical display that simultaneously describes several 
% important features of a data set, such as center, spread, departure from 
% symmetry, and identification of observations that lie unusually far from
% the bulk of the data. In notched Box-plot the notches represent a robust
% estimate of the uncertainty about the medians. If the notches for 2 boxes 
% do not overlap, indicate that the medians of the two groups differ at the
% 5% significance level. However, the comparison may be misleading when the
% notches extend to the end of the box, as they may have been truncated to 
% the end of the box; this may occur with small sample sizes.
%
% Example
%   plotbox({randn(10,1)*5+140, randn(13,1)*8+135});
%   title('Grade 3 heights');
%   set(gca,'xtick',1:2,'xticklabel',{'girls';'boys'});
%   axis([0,3 -inf inf]);
%
% See also plotnorm, plotdensity, plotkde

% Author: Alberto Terruzzi <t-albert@libero.it>
% Version: 1.4
% Created: 6 January 2002

% Version: 1.4.1
% Author: Alberto Pose <apose@alu.itba.edu.ar>
% Updated: 3 September 2006 
% - Replaced deprecated is_nan_or_na(X) with (isnan(X) | isna(X)) 
% (now works with Octave 2.9.7 and foward)
% -revised pab Oct 2007
%  -Translated from Octave to matlab style
%  -renamed from boxplot -> plotbox
% -revised pab June 2008
%  -added options struct
%  - added width to options
%  - reorganized documentation


% Copyright (C) 2002 Alberto Terruzzi
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA


%function s1 = plotbox(data,notched,symbol,vertical,maxwhisker)
function s1 = plotbox(data,varargin)

% assign parameter defaults
if nargin < 1 || nargin > 5
   disp('s = boxplot (data,notch,symbol,vertical,maxwhisker)')
end
options = struct('maxwhisker',1.5,'vertical', 1,'symbol', ['+','o'],'notched', 0,'width','auto');
options = parseoptions(options,varargin{:});
% if nargin < 5 || isempty(maxwhisker), maxwhisker = 1.5; end
% if nargin < 4 || isempty(vertical), vertical = 1; end
% if nargin < 3 || isempty(symbol), symbol = ['+','o']; end
% if nargin < 2 || isempty(notched), notched = 0; end

symbol = options.symbol;
maxwhisker = options.maxwhisker;
vertical = options.vertical;
notched = options.notched;

if length(symbol)==1, symbol(2)=symbol(1); end

if notched==1, notched=0.25; end
a=1-notched;

% figure out how many data sets we have
if iscell(data), 
  nc = length(data);
else
  if isvector(data), data = data(:); end
  nc = size(data,2);
end

% compute statistics
% s will contain
%    1,5    min and max
%    2,3,4  1st, 2nd and 3rd quartile
%    6,7    lower and upper confidence intervals for median
s = zeros(7,nc);
box = zeros(1,nc);
whisker_x = ones(2,1)*[1:nc,1:nc];
whisker_y = zeros(2,2*nc);
outliers_x = [];
outliers_y = [];
outliers2_x = [];
outliers2_y = [];

for i=1:nc
  % Get the next data set from the array or cell array
  if iscell(data)
    col = data{i}(:);
  else
    col = data(:,i);
  end
  % Skip missing data
  col(isnan(col)) = [];
  % Remember the data length
  nd = length(col);
  box(i) = nd;
  if (nd > 1)
    % min,max and quartiles
    s(1:5,i) = [min(col);percentile(col,[0.25; 0.5; 0.75]);max(col)];
    %s(1:5,i) = statistics(col)(1:5);
    % confidence interval for the median
    est = 1.57*(s(4,i)-s(2,i))/sqrt(nd);
    s(6,i) = max([s(3,i)-est, s(2,i)]);
    s(7,i) = min([s(3,i)+est, s(4,i)]);
    % whiskers out to the last point within the desired inter-quartile range
    IQR = (s(4,i)-s(2,i));
    MWR = maxwhisker*IQR;
    whisker_y(:,i) = [min(col(col >= s(2,i)-MWR)); s(2,i)];
    whisker_y(:,nc+i) = [max(col(col <= s(4,i)+MWR)); s(4,i)];
    % outliers beyond 1 and 2 inter-quartile ranges
    outliers = col((col < s(2,i)-MWR & col >= s(2,i)-2*MWR) | (col > s(4,i)+MWR & col <= s(4,i)+2*MWR));
    outliers2 = col(col < s(2,i)-2*MWR | col > s(4,i)+2*MWR);
    outliers_x = [outliers_x; i*ones(size(outliers))];
    outliers_y = [outliers_y; outliers];
    outliers2_x = [outliers2_x; i*ones(size(outliers2))];
    outliers2_y = [outliers2_y; outliers2];
  elseif (nd == 1)
    % all statistics collapse to the value of the point
    s(:,i) = col;
    % single point data sets are plotted as outliers.
    outliers_x = [outliers_x; i];
    outliers_y = [outliers_y; col];
  else
    % no statistics if no points
    s(:,i) = NaN;
  end
end

% Note which boxes don't have enough stats
chop = find(box <= 1);
    
% Draw a box around the quartiles, with width proportional to the number of
% items in the box. Draw notches if desired.
switch options.width
  case 'auto'
    box = box*0.4/max(box);
  case 'fixed'
    box(:) = 0.2;
  otherwise
    box(:) = options.width/2;
end
quartile_x = ones(11,1)*(1:nc) + [-a;-1;-1;1;1;a;1;1;-1;-1;-a]*box;
quartile_y = s([3,7,4,4,7,3,6,2,2,6,3],:);

% Draw a line through the median
median_x = ones(2,1)*(1:nc) + [-a;+a]*box;
median_y = s([3,3],:);

% Chop all boxes which don't have enough stats
quartile_x(:,chop) = [];
quartile_y(:,chop) = [];
whisker_x(:,[chop,chop+nc]) = [];
whisker_y(:,[chop,chop+nc]) = [];
median_x(:,chop) = [];
median_y(:,chop) = [];

% Add caps to the remaining whiskers
cap_size = 0.1; % 0.05
cap_x = whisker_x;
cap_x(1,:) = cap_x(1,:)-cap_size;
cap_x(2,:) = cap_x(2,:)+cap_size;
cap_y = whisker_y([1,1],:);

%#quartile_x,quartile_y
%#whisker_x,whisker_y
%#median_x,median_y
%#cap_x,cap_y

% Do the plot
if vertical
  plot (quartile_x, quartile_y, 'b',whisker_x, whisker_y, 'b',cap_x, cap_y, 'b',median_x, median_y, 'r',...
    outliers_x, outliers_y, [symbol(1),'r'],  outliers2_x, outliers2_y, [symbol(2),'r']);
else
  plot (quartile_y, quartile_x, 'b',	whisker_y, whisker_x, 'b',	cap_y, cap_x, 'b',...
	median_y, median_x, 'r',	outliers_y, outliers_x, [symbol(1),'r'],  outliers2_y, outliers2_x, [symbol(2),'r']);
end
if nargout>0
  s1 = s;
end

end %function
