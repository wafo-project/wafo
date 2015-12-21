%LNKGENPAR Link for x,F and parameters of Generalized Pareto distribution
%  
%   CALL  phati = lnkgenpar(x,logR,phat,i)
%  
%     phati = parameter i as function of x, logR and phat(j) where j ~= i
%     x     = quantile
%     logR  = logarithm of the survival probability
%     
%   LNKGENPAR is a function connecting the quantile (x) and the survival 
%   probability (R) with the fixed distribution parameter, i.e.: 
%     phat(i) = link(x,logR,phat,i), 
%    where logR = log(Prob(X>x;phat)).
%  
%   Example % See proflog
%   
%   See also proflog

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


% Reference 
% Stuart Coles (2004)
% "An introduction to statistical modelling of extreme values".
% Springer series in statistics


function  phati  = lnkgenpar(x,logR,phat,ix)

if numel(phat)<3
  u = 0;
else
  u = phat(3);
end
 
switch ix
  case 1
    error('lnkgenpar(x,logR,phat,i) where i=1 is not implemented!')
  case 2
    % % Reorganizing w.r.t. phat(2) (scale), Eq. 4.13 and 4.14, pp 81 in Coles (2004) gives
    %   link = @(x,logR,phat,ix) -(x-phat(3)).*phat(1)./expm1(phat(1).*logR);
    if phat(1)~=0
      phati =  -(x-u).*phat(1)./expm1(phat(1).*logR);
    else
      phati =  -(x-u)./logR;
    end
  case 3
     if phat(1)~=0
      phati =  x + phat(2).*expm1(phat(1).*logR)./phat(1);
    else
      phati = x+phat(2).*logR;
    end
  otherwise
    error('Index to the fixed parameter is out of bounds')
end