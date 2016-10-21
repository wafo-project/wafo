function [value,bound,inform] = prbnormtndpc(RHO,A,B,D,NDF,abseps,IERC,HNC)
%PRBNORMTNDPC Multivariate normal or T probability with product correlation structure.
%
%  CALL [value,bound,inform] = prbnormtndpc(RHO,A,B,D,NDF,abseps,IERC,HINC)
%
%     RHO    = REAL, array of coefficients defining the correlation
%              coefficient by:
%                correlation(I,J) =  RHO(I)*RHO(J) for J/=I
%              where 
%                1 < RHO(I) < 1
%     A	     = vector of lower integration limits.
%     B	     = vector of upper integration limits.
%	       NOTE: any values greater the 37 in magnitude, 
%                    are considered as infinite values.
%     D      = vector of means (default zeros(size(RHO)))
%     NDF    = Degrees of freedom, NDF<=0 gives normal probabilities (default)
%     ABSEPS = absolute error tolerance. (default 1e-4)
%     IERC   = 1 if strict error control based on fourth derivative
%              0 if intuitive error control based on halving the
%                      intervals (default)
%     HINC   = start interval width of simpson rule (default 0.24)
%
% OUTPUT:
%     VALUE  = estimated value for the integral
%     BOUND  = bound on the error of the approximation
%     INFORM = INTEGER, termination status parameter:
%            0, if normal completion with ERROR < EPS;
%            1, if N > 1000 or N < 1.
%            2, IF  any abs(rho)>=1      
%            4, if  ANY(b(I)<=A(i))
%            5, if number of terms computed exceeds maximum number of
%                  evaluation points
%            6, if fault accurs in normal subroutines
%            7, if subintervals are too narrow or too many
%            8, if bounds exceeds abseps
%
% PRBNORMTNDPC calculates multivariate normal or student T probability
% with product correlation structure for rectangular regions.
% The accuracy is as best around single precision, i.e., about 1e-7.
%   
% Example:
%  rho2 = rand(1,2); 
%  a2   = zeros(1,2);
%  b2   = repmat(inf,1,2);
%  [val2,err2] = prbnormtndpc(rho2,a2,b2);
%  g2 = inline('0.25+asin(x(1)*x(2))/(2*pi)');
%  E2 = g2(rho2)  % exact value
% 
%  rho3 = rand(1,3); 
%  a3   = zeros(1,3);
%  b3   = repmat(inf,1,3);
%  [val3,err] = prbnormtndpc(rho3,a3,b3);  
%  g3 = inline('0.5-sum(sort(acos([x(1)*x(2),x(1)*x(3),x(2)*x(3)])))/(4*pi)');
%  E3 = g3(rho3)   %  Exact value  
%  
% See also  prbnormndpc, prbnormnd, rind
  
% Reference:
%   Charles Dunnett (1989)
%  "Multivariate normal probability integrals with product correlation
%  structure", Applied statistics, Vol 38,No3, (Algorithm AS 251)
    
% The mex-interface and m-file was written by
%     Per Andreas Brodtkorb
%     Norwegian Defence Research Establishment
%     P.O. Box 115
%     N-3191 Horten
%     Norway
%     Email: Per.Brodtkorb@ffi.no
%

% revised pab April 2008
% -renamed from mvnortpcprb to prbnormtndpc

  error(nargchk(3,8,nargin))
  if nargin<4||isempty(D),
    D = zeros(size(RHO));
  end
  if nargin<5||isempty(NDF),
     NDF = 0;
  end
  if nargin<6||isempty(abseps),
     abseps = 1.0e-4;
  end
  if nargin<7||isempty(IERC),
     IERC = 0;
  end
  if nargin<8||isempty(HNC),
     HNC = 0.24;
  end
  % Make sure integration limits are finite
  A = min(max(A,-100),100);
  B = max(min(B,100),-100);
  
  [value,bound,inform] = mvnprdmex(RHO,A,B,D,NDF,abseps,IERC,HNC);
  
