function [len,bin,val] = bincount(x,f)
%BINCOUNT 1-dimensional Bin Count
%
%  CALL: [len,bin, val] = bincount(x,f);
%        [len,bin]      = bincount(x);        
%
%  len = vector with the number of equal values in x, 
%        i.e., len(k) = sum(x==bin(k)).
%  bin = same values as in x, but with no repetitions, 
%        i.e., bin = unique(x).
%  val = vector with the sum of the corresponding values
%        i.e., val(k) = sum(f(x==bin(k))).  
%  x   = vector of function arguments, e.g. an integer index vector.
%  f   = vector of function values, i.e., f(x).
%
% BINCOUNT counts the number of equal values in X, and optionally 
% adds together any elements of F which have duplicate values of X into VAL.
%  
%  Example: 
%  N  = 500; dx = 0.2;
%  f  = rndray(1,N,1);
%  ix = floor(f/dx)+1;
%  [len,bin] = bincount(ix); 
%  plot((bin-.5)*dx,len/N/dx,'.'); % 1D probability density plot
%  bar((bin-.5)*dx,len/N/dx);      % 1D probability density plot
%  bar((bin-.5)*dx,len);           % 1D Histogram
%
%  close all;
%
% See also  sparse, histc.

%Tested on: Matlab 5.3
%History:
% revised pab jan2007
% -replace call to csort with sort since it is faster.
% revised pab Dec2003
% renamed from binc to bincount  
% by pab 14.08.2001

% check number of input arguments
error(nargchk(1, 2, nargin));
if isempty(x)
  len = 0;
  bin = [];
  val = 0;
  return
end
isiz = size(x);
ldim = isiz > 1;
if sum(ldim) > 1
  error('Input must be a vector.');
end

% make sure input is a column vector
x = x(:);

[x, ind] = sort(x);

% Find indices to unique values
i = [ find(diff(x)) ; length(x) ];

if nargout>1,   
  bin = x(i); % bin = unique(x);
end 
i = [ 0 ; i ];

len = diff(i);
   
% Make sure that the output is a vector in the same dimension as input
osiz       = isiz;
osiz(ldim) = length(len);
len        = reshape(len, osiz);
   
val = [];
if (nargin>1 && nargout>2),
  if any(isiz~=size(f)),
    error('The size of x and f must be equal!'),
  end
  f   = f(:);
  f   = [ 0; cumsum(f(ind))];
  val = diff(f(i+1));
  val = reshape(val, osiz);
end







