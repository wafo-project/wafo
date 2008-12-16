function phat = plotgumb(x)
%PLOTGUMB Plot data on Gumbel distribution paper.
%
% CALL:  phat = plotgumb(X)
%
%       phat = [a b] Parameters (see prbgumb) estimated from the plot by
%              least squares method 
%          X = data vector or matrix
%
% Example:
%   R=rndgumb(2,0,1,100);
%   phat=plotgumb(R),shg
%
% See also pdfgumb, prbgumb, rndgumb, fitgumb, momgumb

% Reference: 
%  Johnson  N.L., Kotz S. and Balakrishnan, N. (1994)
%  Continuous Univariate Distributions, Volume 2. Wiley. 


% rewritten ms 20.06.2000

F=edf(x,'wdata',true);
x = F.args(:);
F1 = F.data;
plot(x,-log(-log(F1)),'b.','markersize',12);
U=[ones(size(x)) x];
c=U\(-log(-log(F1)));
a=1/c(2);
b=-c(1)*a;
hold on
plot(x,U*c,'r--')
hold off
title('Gumbel Probability Plot')
xlabel('X')
ylabel('-log(-log(F))')
wafostamp;
if nargout > 0,
  phat=[a,b];
end
