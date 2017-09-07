function h1=tallibing(x,y,n,size,color)
% TALLIBING  Display numbers on field-plot
%
% CALL h=tallibing(x,y,n,size,color)
%
% x,y    = position matrices
% n      = the corresponding matrix of the values to be written
%          (non-integers are rounded)
% size   = font size (optional) (default=8)
% color  = color of text (optional) (default='white')
% h      = column-vector of handles to TEXT objects
%
% TALLIBING writes the numbers in a matrix as text at the positions 
% given by the x and y coordinate matrices.
%  When plotting binned results, the number of datapoints in each
%          bin can be written on the bins in the plot.
%
% EXAMPLE: 
%  [x,y,z] = peaks(20); 
%  epcolor(x,y,z); 
%  tallibing(x,y,z);
%   % pcolor(x,y,z); shading interp; 
%
%  close all;
%  
% See also TEXT TALLIBING3

%Time-stamp:<Last updated on 01/01/07 at 12:27:30 by even@gfi.uib.no>
%File:<d:/home/matlab/tallibing.m>

%error(nargchk(3,5,nargin));
narginchk(3,5)
if nargin < 5 || isempty(color), color='w'; end
if nargin < 4 || isempty(size),  size=8;    end

if isvector(x) || isvector(y)
  [x,y]=meshgrid(x,y);
end

x=x(:); y=y(:); n=n(:);
n=round(n);

oldtall=findobj(get(gca,'children'),'Tag','tallibing');
if any(oldtall), delete(oldtall); end
h = cell(1,0);
for i=1:length(y)
  if n(i)
    h{end+1}=text(x(i),y(i),int2str(n(i)),...
           'HorizontalAlignment','center','FontWeight','demi',...
           'FontSize',size,'Color',color,...
	   'Tag','tallibing');
  end
end

if nargout>0
  h1 = h;
end
