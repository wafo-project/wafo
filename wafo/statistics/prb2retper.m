function retper = prb2retper(prb,npu)
%PRB2RETPER Return period from Probability of exceedance.  
%
% CALL retper = retper2prb(prb,npu)
%
%   retper  = the return period
%   prb     = probability of not exceeding level
%   npu     = the mean Number of events Per Unit time (eg Year)
%
% Example
%
% See also retper2prb

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


  if (any(npu) <=0)
    error('NPU must be positive!')
  end
  if (any(prb<0 | 1 <= prb))
    error('Probability must be the range [0, 1)!')
  end
  
  
  retper = 1 ./ (npu * (1 - prb));
  
