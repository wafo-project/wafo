function prb = retper2prb(retper,npu)
%RETPER2PRB Probability of exceedance from return period. 
%
% CALL prb = retper2prb(retper,npu)
%
%   prb     =  probability of exceeding level
%   retper  = the return period
%   npu     = the mean Number of events Per unit time (eg Year)

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
  if (any(retper < 1/npu))
    error('return period incompatible with NPU!')
  end
  prb =  1./(npu .* retper);
  
