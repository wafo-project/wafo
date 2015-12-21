function h1 = labelfig(f,ax)
%WDATA/LABELFIG Labels the figure with title, caption, x-,y- and/or -zlabel
% 
% CALL labelfig(f)
%
% f = wdata object
%
% 
% See also wafostamp

% Tested on: matlab 7
% History:
% By pab Jan 2007



h2 = [];
if nargin>1
  oldAx = gca;
  % make axes ax current
  axes(ax)
end

if ~isempty(f.labels)
  hf   = {@xlabel, @ylabel, @zlabel};
  Nlab = length(f.labels);
  h2   = zeros(1,Nlab);
  for ix = 1:Nlab
    h2(ix) = feval(hf{ix},f.labels{ix});
  end
end
if ~isempty(f.title)
  h2 = [title(f.title),h2];
end
if nargout>0
  h1 = h2;
end

caption = '';
if ~isempty(f.caption)
  caption = f.caption;
end
wafostamp(caption)

if exist('oldAx','var')
  axes(oldAx)
end