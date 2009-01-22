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
