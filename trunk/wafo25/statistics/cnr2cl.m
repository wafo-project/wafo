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
