function r = ranktrf(x)
%RANKTRF  Rank transformation of data material.
%
%	CALL  r = ranktrf(x)
%
%  Replaces the data vector  x  and gives rank numbers.
%	 In case  x  is a matrix every column is rank transformed.
%  
%  Ties in the data are resolved by using average rank value
%
% Example
%  a = rand(5,1);
%  r = ranktrf(a(:));
%  [a,r]
%  a2 = a;
%  a2(1) = a(2);
%  r2 = ranktrf(a2(:));
%  [a2,r2]
% 
% See also spearman, testmean1r, testmean2r

%      Copyright (c) Anders Holtsberg, 1998

tr = 0;
if size(x,1) == 1
   x = x';
   tr = 1;
end

%n = size(x,1);
%v = (1:n)';

[dummy, i]  = sort(x);
[dummy, r1] = sort(i);
[dummy, i]  = sort(flipud(x));
[dummy, r2] = sort(i);
r = (r1+flipud(r2))/2;

if tr
   r = r';
end
