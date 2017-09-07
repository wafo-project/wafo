function s = poly2hstr( coefs, var )
%POLY2HSTR Polynomial as a Horner represented string.
%
%   S = POLY2HSTR( P ) returns a string S consisting of the polynomial
%   coefficients in the vector P multiplied by powers of the variable
%   'x'.
%
%   S = POLY2HSTR( P, 's' ) uses 's' as the variable rather than 'x'.
%
%  Example: 
%  assert(poly2hstr( [1 1 2], 's' ), '(s + 1)*s + 2');
%
% See also poly2str

% History
% by pab 16.02.2001
% based on poly2str by
%   Author:      Peter J. Acklam
%   Time-stamp:  1998-06-22 20:37:43
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

%error( nargchk( 1, 2, nargin ) );
narginchk(1,2)
if nargin == 1
   var = 'x';
end
coefs = polytrim(coefs);
order = length(coefs)-1;        % Order of polynomial.
s = '';                         % Initialize output string.

ix=1;
for expon = order:-1:0
   coef = coefs( order-expon+1 );
   % There is no point in adding a zero term (except if it's the only
   % term, but we'll take care of that later).
   if coef ~= 0
     % Append exponent if necessary.
     if ix>1
        exponstr = sprintf( '%.0f', ix);
        s = [ s '^' exponstr ];
        ix=1;
     end

     % Is it the first term?
     isfirst = isempty(s);
     
     % We need the coefficient only if it is different from 1 or -1 or
     % when it is the constant term.
     needcoef = (( abs(coef) ~= 1 ) | ( expon == 0 ) & isfirst) | ~isfirst;
     
     % We need the variable except in the constant term.
     needvar = ( expon ~= 0 );
          
     % Add sign, but we don't need a leading plus-sign.
     if isfirst
       if coef < 0
	 s = [ '-' ];        % Unary minus.
       end
     else
       if coef < 0
	 s = [ s ' - ' ];    % Binary minus (subtraction).
       else
	 s = [ s ' + ' ];    % Binary plus (addition).
       end
     end
     
     % Append the coefficient if it is different from one or when it is
     % the constant term.
     if needcoef
       coefstr = sprintf( '%.20g', abs(coef) );
       %coefstr = sprintf( '%.15g', abs(coef) );
       s = [ s coefstr ];
     end
     
     
     % Append variable if necessary.
     if needvar
       % Append a multiplication sign if necessary.
       if needcoef 
	 if ~isfirst, 
	   s = ['(', s ,')' ];
	 end
	 s = [ s ,'*' ];
       end
       s = [ s var ];
     end
   else
      ix = ix+1;
   end
 end

% Now treat the special case where the polynomial is zero.
if isempty(s)
   s = '0';
end
