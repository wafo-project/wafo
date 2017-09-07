function s = vectoriz(s)
%VECTORIZ Vectorize expression or inline function object.
%
% CALL: v = vectoriz(s);
% 
%  s = String expression or inline object
%  v = Vectorized expression or inline object.
%
% VECTORIZ inserts a '.' before any '^', '*' or '/' in S  
% not already proceeded by a '.'.
%
% Example:
% assert(vectoriz('x^2 + 1/x + x*y'), 'x.^2 + 1./x + x.*y')
%
% See also: inline/formula,


%History
% by pab 25.01.2001

%error(nargchk(1,1,nargin));
narginchk(1,1)
switch class(s)
  case 'inline',  v = s.expr; %formula(s); 
    islin=1;
  case 'char',    v = s;    islin=0;
  otherwise, 
    error('VECTORIZ: Input must be string expression or inline object');
end

ind = find((v=='^') | (v=='*') | (v=='/'));
if any(ind),
  ind(v(max(ind-1,1))=='.')=[]; % remove indices to parts already vectorized
  v = insert(v,'.',ind);
  if islin==1,
    s.expr = v;
  else
    s = v;
  end
end
return