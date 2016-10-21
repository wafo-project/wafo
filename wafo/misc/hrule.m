function [bp,wf] = hrule(n,method)
%HRULE  Computes nodes and weights for generalized Hermite quadrature.
% 
%  CALL: [bp,wf] = hrule(n,method)
%
%  bp = base points
%  wf = weight factors
%  n  = number of base points (integrates a (2n-1)th order
%       polynomial exactly)
% method = Integer defining the method used to obtain the nodes and weights
%          0 Newton Raphson to find zeros of the Hermite polynomial (Fast)
%            (Default)
%          1 Eigenvalue of the jacobi matrix (Slow)
% 
%
%  The Hermite Quadrature integrates an integral of the form
%       inf                          n
%      Int (exp(-x.^2) F(x)) dx  =  Sum  wf(j)*F( bp(j) )
%     -inf                           j=1    
%
% Example
%  [x,w]= hrule(10)
%  sum(x.*w)
% 
%  See also qrule, gaussq

% Reference 
%   method 0: see ref [4]
%   method 1:  Adapted from Netlib routine gaussq.f see ref [2,3]
%
% [2]  Golub, G. H. and Welsch, J. H. (1969)
% 'Calculation of Gaussian Quadrature Rules'
%  Mathematics of Computation, vol 23,page 221-230,
%
% [3]. Stroud and Secrest (1966), 'gaussian quadrature formulas', 
%      prentice-hall, Englewood cliffs, n.j.
% 
% 
% [4]William H. Press, Saul Teukolsky, 
% William T. Wetterling and Brian P. Flannery (1997)
% "Numerical recipes in Fortran 77", Vol. 1, pp 145
% "Numerical recipes in Fortran 90"

% Adopted to matlab by pab jan 2006

if nargin<2 || isempty(method)
  method =0;
end

switch method
  case 1,
    ii = (1:(n-1))';
    muzero = sqrt(pi);
    a      = zeros(n,1);
    b      = sqrt(ii/2);
    [bp,wf] = gausscof(a,b,muzero);
  otherwise

    MAXIT=10;
    releps = 3e-14;
    C=[9.084064e-01;...
      5.214976e-02;...
      2.579930e-03;...
      3.986126e-03];
    %PIM4=0.7511255444649425_dp
    PIM4 = pi^(-1/4);
    
    % The roots are symmetric about the origin, so we have to
    % find only half of them.
    m=fix((n+1)/2) ;
    
    % Initial approximations to the roots go into z.
    anu = 2.0*n+1;
    rhs = (3:4:4*m-1)*pi/anu;
    r3  = rhs.^(1/3);
    r2  = r3.^2;
    theta = r3.*(C(1)+r2.*(C(2)+r2.*(C(3)+r2.*C(4))));
    z = sqrt(anu).*cos(theta);
    
    L = zeros(3,length(z));
    k0  = 1;
    kp1 = 2;
    for its=1:MAXIT
      %Newton’s method carried out simultaneously on the roots.
      L(k0,:)  = 0;
      L(kp1,:) = PIM4;
      
      for j=1:n
        %Loop up the recurrence relation to get the Hermite
        %polynomials evaluated at z.
        km1 = k0;
        k0  = kp1;
        kp1 = mod(kp1,3)+1;
        
        L(kp1,:) =z.*sqrt(2/j).*L(k0,:)-sqrt((j-1)/j).*L(km1,:);
        
      end
      %L now contains the desired Hermite polynomials.
      % We next compute pp, the derivatives,
      % by the relation (4.5.21) using p2, the polynomials
      % of one lower order.
      
      pp = sqrt(2*n)*L(k0,:);
      dz = L(kp1,:)./pp;
      
      z=z-dz; % Newton’s formula.
      unfinished=find(abs(dz) > releps);
      
      if (not( any(unfinished)))
        break;
      end
    end
    if any(unfinished)
      warning('WAFO:HRULE','too many iterations in hrule');
    end
    %bp = [z,-z];
    bp(1:m)        = z;      % Store the root
    bp(n:-1:n-m+1) = -z;     % and its symmetric counterpart.
    wf(1:m)        = 2./pp.^2; % Compute the weight
    wf(n:-1:n-m+1) = wf(1:m);  % and its symmetric counterpart.
end