%LNKRAY Link for x,F and parameters of Rayleigh distribution
%  
%   CALL  phati = lnkray(x,logR,phat,i)
%  
%     phati = parameter i as function of x, logR and phat(j) where j ~= i
%     x     = quantile
%     logR  = logarithm of the survival probability
%     
%   LNKRAY is a function connecting the quantile (x) and the survival 
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




function  phati  = lnkray(x,logR,phat,ix)

phati =  x./ sqrt(-2*logR);
