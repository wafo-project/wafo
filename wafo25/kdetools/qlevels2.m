function [q, p]=qlevels2(r,p,method)
%QLEVELS2 Calculates quantile levels which encloses P% of data
%
%  CALL: [ql PL] = qlevels2(data,PL,method);
%
%    ql   = the discrete quantile levels, size Np x D
%    data = data matrix, size N x D (D = # of dimensions)
%    PL   = percent level vector, length Np (default [10:20:90 95 99 99.9])
%  method = 1 Interpolation so that F(X_(k)) == (k-0.5)/n. (default)
% 	    2 Interpolation so that F(X_(k)) == k/(n+1).
% 	    3 Based on the empirical distribution.
%
% QLEVELS2 sort the columns of data in ascending order and find the  
%          quantile levels for each column which encloses  P% of the data.  
%  
% Examples : % Finding quantile levels enclosing P% of data: 
%   xs  = rndnorm(0,1,100000,1);
%   qls = qlevels2(pdfnorm(xs),[10:20:90 95 99 99.9]);
%            % compared with the exact values
%   ql  = pdfnorm(invnorm((100-[10:20:90 95 99 99.9])/200));
%
% % Finding the median of xs:
%   ql  = qlevels2(xs,50); 
%
% See also  qlevels

%tested on: matlab 5.3
% History:
% revised pab 
%  -fixed a bug for D>=3
%  -added different methods from wquantile
%  -updated example  
% by pab 25.09.1999
% 

[n, d]=size(r);
if nargin<2||isempty(p)
  p = [10:20:90 95 99 99.9];
elseif any(p<0 | 100<p),
  error('PL must satisfy 0 <= PL <= 100')
end
if nargin<3||isempty(method), method=1; end

if (n==1) && (d>1)
  r=r(:);
  n=d;
  d=1;
end
if d>1
  if min(size(p)) > 1 
    error('Not both matrix r and matrix p input')
  end
  q = zeros(length(p),d);
else
  q = zeros(size(p));
end
p = 1-p(:)/100;
x = sort(r);


if method == 3
   qq1 = x(ceil(max(1,p*n)),:); 
   qq2 = x(floor(min(p*n+1,n)),:);
   qq = (qq1+qq2)/2;
else                         
   x = [x(1,:); x; x(n,:)];
   if method == 2
      % This method is from Hjort's "Computer
      % intensive statistical methods" page 102
      i = p*(n+1)+1;
   else % Metod 1
      i = p*n+1.5;
   end
   iu = ceil(i);
   il = floor(i);
   d1 = (i-il)*ones(1,d);
   qq = x(il,:).*(1-d1)+x(iu,:).*d1;
end

q(:) = qq;

return


