function [LL,pcov,H] = likgenpar(phat,data)
%LIKGENPAR Log likelihood function for GPD
%
% CALL  [L,pcov] = likgenpar(phat,data)
%
% Example
%
%
% See also fitgenpar

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


n = numel(data);
k = phat(1);
s = phat(2);
m = 0;
np = length(phat);
if np>2
  m = phat(3);
end
xn = (data-m)./s;
if min(xn)<0 || (k>0 && k*max(xn)>=1)
  LL = Inf;
  if nargout > 1
      H = [NaN NaN; NaN NaN];
    pcov = H;
  end
elseif abs(k) > eps
  kxn = k.*xn;
  sumlog1mkxn = sum(log1p(-kxn));
  LL = n*log(s) + (1-1/k).*sumlog1mkxn;
  if nargout > 1
    r = xn./(1-kxn);
    sumr = sum(r);
    sumr2 = sum(r.^2);
    H11 = 2*sumlog1mkxn./k^3 + 2*sumr./k^2 + (1-1/k).*sumr2;
    H12 = -(sumr + (k-1).*sumr2)./s;
    H22 = (n + 2*(k-1).*sumr + k*(k-1).*sumr2)./s^2;
    H = [H11,H12;H12,H22];
    % Invert observed information number
    pcov = -[H11 H12; H12 H22]\eye(2);
  end
else
  % k == 0
  sumx = sum(xn);
  LL = n*log(s) + sumx;
  if nargout > 1
    sumx2 = sum(xn.^2);
    H11 = -(2/3)*sum(xn.^3) + sumx2;
    H12 = -(sumx - sumx2)./s;
    H22 = (n - 2*sumx)./s^2;
    H = [H11,H12;H12,H22]; % Hessian matrix
    pcov = -[H11,H12;H12,H22]\eye(2);
  end
end