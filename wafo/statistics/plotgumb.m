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

%
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU Lesser General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.


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
