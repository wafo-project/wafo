function I1 = cnr2cl(I)
%CNR2CL Column Number to Column Label.
%
% CALL:  I1 = cnr2cl(I)
%
% I1 = character array of column labels.
% I  = matrix of column numbers
%
% CNR2CL transforms a column number into a column label, i.e., 
% convert 'A'-'Z' to 1-25, 'a'-'z' to 26-50 and ' ' to 0.
%
% CNR2CL is useful in conjuction with SUDG and CDR
% 
% See also  cl2cnr, sudg, cdr

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

if isnumeric(I)
  if any(I(:)>50),
    warning('WAFO:CNR2CL','Column numbers must be less than 51!')
    I1 = I;
  else
    
    % characters ' ','A' - 'Z' 'a'-'z'
    str1=[' ',char(65:90) char(97:122)];
    I1  = str1(abs(I)+1);
    sgn         = sign(I);
    sgn(sgn==0) = 1; % define sign(0) to 1.
    sgn         = prod(sgn,2);
    k = find(sgn<0);
    if any(k),  % add '-' sign to strings
      sgn(k) = 2;             % Define negative values to 2.
      str0 = ' -'; % '+' and '-'
      I1   = [str0(sgn).' I1];
    end
  end
else
  I1 = I;
end

return
