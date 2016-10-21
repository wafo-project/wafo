function phat = plotexp(x)
%PLOTEXP Plot data on Exponential distribution paper
%
% CALL:  phat = plotexp(X)
%
%       phat = [m] Parameter (see prbexp) estimated from 
%              the plot by least squares method
%          X = data vector or matrix
%
% Example:
%   R=rndexp(2,1,100);
%   phat=plotexp(R),shg
%
% See also  cdfexp, plotweib

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


x = x(:);
F=edf(x,'wdata',true');
x  = F.args;
F1 = F.data;
plot(x,-log1p(-F1),'b.','markersize',12);

m = mean(x);
hold on
plot(x,1/m*x,'r--')
hold off
title(['Exponential Probability Plot, m=' num2str(m)])
xlabel('x')
ylabel('-log(1-F)')
if nargout > 0,
  phat= m;
end
wafostamp;
