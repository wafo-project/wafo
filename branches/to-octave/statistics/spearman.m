function rho = spearman(x,y)
%SPEARMAN Spearman's rank correlation coefficient.
%
%	  rho = spearman(x,y)
%
%	  This is the correlation coefficient between rank
%	  transformed data. 

if nargin < 2
   rho = corrcoef(ranktrf(x));
else
   rho = corrcoef(ranktrf(x),ranktrf(y));
end