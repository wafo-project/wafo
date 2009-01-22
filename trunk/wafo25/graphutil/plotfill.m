function h = plotfill(x,y,color,varargin)
%PLOTFILL
% 
% Example
%  x = linspace(0,10).';
%  y = [sin(x)+1,sin(x)-1];
%  plotfill(x,y,'g','facealpha',0.1)
% 
% See also fill
error(nargchk(2,inf,nargin))
[n,m] = size(y);
if (n==2 || m==2) && (n>=2 && m>=2)
  if n==2 && m>2;
    y =y.';
    [n,m] = size(y);
  end
else
  error('Y must have size 2 x M or N x 2!')
end
nx = numel(x);
if nx~=n
  error('The length of X must equal the length of Y!')
end
if nargin<3 || isempty(color)
  color = 'r';
end
%x = linspace(0,10).';
x2 = x(end:-1:1);
%y = sin(x);
y1 = y(:,1);
y2 = y(n:-1:1,2);

h1 = fill([x;x2],[y1;y2],color,'facealpha',1,'FaceVertexAlphaData',0.5,'edgealpha',0,varargin{:});

if nargout>0
  h = h1;
end