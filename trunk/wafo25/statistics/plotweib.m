function phat = plotweib(x)
%PLOTWEIB Plot data on a Weibull distribution paper
%
% CALL:  phat = plotweib(X)
%
%       phat = [a c] Parameters (see prbweib) estimated from 
%              the plot by least squares method
%          X = data vector or matrix
%
% Example:
%   R=rndweib(2,2,0,1,100);
%   phat=plotweib(R),shg
%
% See also  cdfweib, invweib

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



% Reference: Cohen & Whittle, (1988) "Parameter Estimation in Reliability
% and Life Span Models", p. 25 ff, Marcel Dekker.

% Revised jr 28.08.2000: line 23 added.
% rewritten ms 20.06.2000

x = x(:);
F=edf(x,'wdata',true);
x = F.args(:);
F1 = F.data;
plot(log(x),log(-log1p(-F1)),'b.','markersize',12);
U=[ones(size(x)) log(x)];
b=U\log(-log1p(-F1));
c=b(2);
a=exp(-b(1)/c);
hold on
plot(log(x),U*b,'r--')
hold off
title('Weibull Probability Plot')
xlabel('log(X)')
ylabel('log(-log(1-F))')
if nargout > 0,
  phat=[a,c];
end

wafostamp;
