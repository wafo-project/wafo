function h1 = tallibing3(x,y,z,n,size,color)
% TALLIBING3 Display numbers in 3d plot
% Writes the numbers in an array as text at the 3D-positions given by the
% x, y and z coordinate-arrays.
%
% tallibing3(x,y,z,n,size,color)
%
% x,y,z  = position 3D-arrays
% n      = the corresponding array of the values to be written
%          (non-integers are rounded)
% size   = font size (optional) (default=8)
% color  = color of text (optional) (default='white')
% h      = column-vector of handles to TEXT objects
%
% USAGE:   When plotting binned results, the number of datapoints in each
%          bin can be written in the bins in the plot.
% Example
%  [x,y,z] = peaks(20); 
%  surf(x,y,z); tallibing3(x,y,z,z);
%
%  close all;
%
% See also TALLIBING TEXT

%Time-stamp:<Last updated on 00/07/11 at 13:36:11 by even@gfi.uib.no>
%File:<d:/home/matlab/tallibing3.m>

%error(nargchk(4,6,nargin));
narginchk(4,6)
if nargin < 6 || isempty(color), color='w'; end
if nargin < 5 || isempty(size),  size=8;    end

x=x(:); y=y(:); z=z(:); n=n(:);
n=round(n);
h = cell(1,0);
for i=1:length(y)
  if n(i)
    h{end+1}=text(x(i),y(i),z(i),int2str(n(i)),...
           'HorizontalAlignment','center','FontWeight','demi',...
           'FontSize',size,'Color',color);
  end
end
if nargout>0
    h1 = h;
end

