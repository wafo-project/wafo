function data = get(Hs,property)
%GET Access data stored in a SPECDATA object
% 
% CALL  data = get(Hs,property)

% Revised pab 01.11.2003
% vectorized code to speed up things  
  

%error(nargchk(2,2,nargin))
narginchk(2,2)
index.type = '.';
index.subs = property;
data = subsref(Hs,index);
%data = Hs.(property);
