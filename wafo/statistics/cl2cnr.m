function I1 = cl2cnr(I)
%CL2CNR Column Label to Column Number
%
% CALL:  I1 = cl2cnr(I)
%
% I1 = matrix of column numbers
% I  = character array of column labels.
%
% CL2CNR transforms a column label into a column number, i.e., 
% convert 'A'-'Z' to 1-25, 'a'-'z' to 26-50 and ' ' to 0.
% 
% CL2CNR is useful in conjuction with SUDG and CDR
%
% See also  cnr2cl, sudg, cdr

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


% Tested on: Matlab 5.3
% History:
% By Per A. Brodtkorb 16.03.2001

if ischar(I)
  sgn = ones(size(I));
  sgn(I=='-') = -1;
  sgn = prod(sgn,2);
  
  I1 = double(I)-64;     % Convert A-Z to 1-25
  I1(I1<0)=0;            % Convert ' ' to 0
  k = find(I>=97);
  if any(k),             % Convert a-z to 26-50
    I1(k) = I(k)-96+26;
  end
  I1 = sort(I1,2);
  
  if any(I1(:)>50), warning('WAFO:CL2CNR','Illegal column label!'), end
  I1(:,end) = I1(:,end).*sgn;
else
  I1 = I;
end

% Remove starting zeros
[ix,iy] = find(I1~=0); 
iy  = min(iy);
if ~isempty(iy) && (iy>1),  I1 = I1(:,iy:end);  end
  
return
