function h=epcolor(x,y,data)
%EPCOLOR  Pseudocolor (checkerboard) plot with mid-bin positioning.
%
% CALL:   h = epcolor(x,y,data)
% 
% h    = handle to surface object
% [x,y]= the axes corresponding to the data-positions. Vectors or
%        matrices. If omitted, giving only data-matrix as inargument, the
%        matrix-indices are used as axes.
% data = data-matrix
%
% EPCOLOR make a checkerboard plot where the data-point-positions are in
% the middle of the bins instead of in the corners, and the last column
% and row of data are used.
%
% Since the x- and ydata of surface object is shifted to fit PCOLOR, the
% centered x and y data is stored in object's UserData as 1x2 cell array.
%
% Example: 
% [x,y,z]=peaks(20); 
% epcolor(x,y,z);
%
% close all;
%
% See also PCOLOR

%Time-stamp: <2000-06-29 18:43:42 even>
%File:</home/even/matlab/evenmat/epcolor.m>
% Revised pab: M
  
  
%error(nargchk(1,3,nargin));
narginchk(1,3)
if nargin==1
  data=x;
  [M,N]=size(data);
  x=1:N;
  y=1:M;
elseif nargin==3
  [M,N]=size(data);
  if min(size(x))~=1, x=x(1,:);  end
  if min(size(y))~=1, y=y(:,1)'; end
else
   error('EPCOLOR takes 3 or 1 inarguments! (x,y,data) or (data)')
end

% getting the right axes (askew [-dx,-dy] + add one at the ends)
[tmp,xx]=buildgrid(x);
[tmp,yy]=buildgrid(y);

% getting datamatrix right (add NaNs in extra row and column)
data(M+1,:)=nan;
data(:,N+1)=nan;

h=pcolor(xx,yy,data);

set(h,'userdata',{x,y});
addtag('epcolor',h);

end

function [X,XG]=buildgrid(X)
% BUILDGRID     builds two vectors describing a 1D-bin structure
%
% [X,XG] = buildgrid(X,x)
%
% Input:
% x  = vector or matrix of the positions of the data that are to be
%      binned
% X  = bin-specifications interpreted as follows:
% 
%      single real number       =  the number of bins (N)
%      single imaginary number  =  the bin-size (dx)
%      real valued vector       =  the bins' mid-point positions
%      imaginary valued vector  =  the limits of the bins
%
%      The first two options result in an evenly spaced grid over the range
%      of x (here x needs to be given).  The last two options gives the
%      bin-structure explicitly either by its mid-point positions or by its
%      limits, and the values of the other is distributed between them. This
%      is useful for making inhomogeneous bin-structures like
% 
%      |.| . |  .  |   .   |    .    |   (dots=positions; lines=limits)
%
%      The imaginary input is simply used as an option-flag. Just multiply
%      your input value by the imaginary unit i.  If X is given as a matrix,
%      the matrix must be invariate in one dimension. Then the other
%      dimension, along which the values change, is used as bin-specs.
%
%      Matrix input is only possible with regular-position-matrices
%      (matrices with repeated rows or columns).
%
% X  = length N vector of the bins' mid-point positions
% XG = length N+1 vector of the limits of the bins
%
% See also BIN1D BIN2D BIN3D MAT2VEC

%Time-stamp:<Last updated on 06/04/24 at 13:19:50 by even@nersc.no>
%File:</home/even/matlab/evenmat/buildgrid.m>

%error(nargchk(1,2,nargin));		% INPUT-TESTS:
narginchk(1,2)
if isempty(X)
  X=20;
end
%x=x(:);				% secure that x is a vector
X=mat2vec(X);			% secure that X is a vector

sx = size(X);
if all(sx==1)
 
  error(['Single valued bin-specification needs to be ',...
    'followed by a vector or regular-grid-matrix ',...
    'of data-point positions!']);
  
else
  
  if isreal(X),		
    XG=findbins(X);
  else
    XG=imag(X);
    X=findpoints(XG);
  end
end
end

%---------------------------------------------------------------------
function X=findpoints(XG)
% X = points half way between all values in XG
X=uplus(XG(1:end-1))+uplus(diff(XG)/2); % preserves shape
end
%---------------------------------------------------------------------
function XG=findbins(X)
% XG = points half way between all values of X _and_ outside the
% endpoints. The outer limits have same distance from X's endpoints as
% the limits just inside.
%dim=isvec(X);			% preserves shape
[N,M] = size(X);
dim = find([N,M]~=1);
dx = diff(X)/2; 
dx = cat(dim,dx,dx(end));%[ans ans(end)];
XG=uplus(X)+uplus(dx);
XG=cat(dim,X(1)-(XG(1)-X(1)),XG);%[X(1)-(XG(1)-X(1)) XG];

end
function x = mat2vec(X)

sx = size(X);
isMatrix = sum(any(sx>1))>1;

if ~isMatrix
  x = X;
  return
end

dx = any(diff(X(1,:)));
dy = any(diff(X(:,1)));

if dx && ~dy
  x = X(1,:);
elseif ~dx && dy
  x = X(:,1);
else
  error('Input matrix is not uniform in any dimension!')
end
end