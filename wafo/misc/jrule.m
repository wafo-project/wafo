function [bp,wf]=jrule(n,alpha,beta,method)
%JRULE  Computes nodes and weights for Gauss-Jacobi quadrature.
% 
%  CALL: [bp,wf] = jrule(n,alpha,beta,method)
%
%  bp = base points
%  wf = weight factors
%  n  = number of base points (integrates a (2n-1)th order
%       polynomial exactly)
% alpha,beta = 
% method = Defines the method used to obtain the nodes and weights
%          0 Newton Raphson to find zeros of the Jacobi polynomial
%          1 Eigenvalue of the jacobi matrix
% 
%
%  The Gauss-Jacobi Quadrature integrates an integral of the form
%       b                                       n
%      Int (x-a)^alpha*(b-x)^beta F(x)) dx  =  Sum  wf(j)*F( bp(j) )
%       a                                      j=1    
%
% Example
%  [x,w]= jrule(10,0,0);
%  
%  assert(sum(x.*w), 2.77555756156289e-016, eps)
%  assert(x, [ 0.973906528517172, 0.865063366688985, 0.679409568299024,...
%              0.433395394129247, 0.148874338981631, -0.148874338981631,...
%             -0.433395394129247, -0.679409568299024, -0.865063366688985,...  
%             -0.973906528517172], 1e-10)
%  assert(w, [0.0666713443086734, 0.1494513491505790, 0.2190863625159803,...
%             0.2692667193099939, 0.2955242247147502, 0.2955242247147502,...
%             0.2692667193099939, 0.2190863625159799, 0.1494513491505790,...   
%             0.0666713443086734], 1e-10)
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

% By pab jan2006
if nargin<4 || isempty(method)
  method = 0;
end

  switch method
    case 1,  %jacobi
      a=zeros(n,1);
      b=zeros(n-1,1);
      ab = alpha + beta;
      abi = 2 + ab;
      muzero = 2^(ab + 1) * gamma(alpha + 1) * gamma(beta + 1) / gamma(abi);
      a(1) = (beta - alpha)/abi;
      b(1) = sqrt(4*(1 + alpha)*(1 + beta)/((abi + 1)*abi^2));
      a2b2 = beta^2 - alpha^2;
      
      i = (2:n-1)';
      abi = 2*i + ab;
      a(i) = a2b2./((abi - 2).*abi);
      a(n) =a2b2./((2*n - 2+ab).*(2*n+ab));
      b(i) = sqrt (4*i.*(i + alpha).*(i + beta).*(i + ab)./((abi.^2 - 1).*abi.^2));
      
      [bp,wf] = gausscof(a,b,muzero);
    otherwise
      MAXIT  = 10;
      releps = 3e-14;
  
      
      % Initial approximations to the roots go into z.
      alfbet = alpha+beta;
    
    
      z = cos( pi*((1:n) -0.25+0.5*alpha)/( n +0.5 *(alfbet+1) ));
      
      L = zeros(3,length(z));
      k0  = 1;
      kp1 = 2;
      for its=1:MAXIT
        %Newton�s method carried out simultaneously on the roots.
        tmp = 2 + alfbet;
        L(k0,:)  = 1;
        L(kp1,:) = (alpha-beta+tmp*z)/2;
        
        for j=2:n
          %Loop up the recurrence relation to get the Jacobi
          %polynomials evaluated at z.
          km1 = k0;
          k0 = kp1;
          kp1 = mod(kp1,3)+1;
          
          a = 2*j*(j+alfbet).*tmp;
          tmp = tmp +2;
          c = 2*(j-1+alpha)*(j-1+beta)*tmp;
          
          b = (tmp-1).*(alpha.^2-beta.^2+tmp.*(tmp-2).*z);
          
          L(kp1,:) =(b.*L(k0,:)-c*L(km1,:))./a;
        end
        %L now contains the desired Jacobi polynomials.
        %We next compute pp, the derivatives with a standard
        % relation involving the polynomials of one lower order.
        
        pp=(n*(alpha-beta-tmp.*z).*L(kp1,:)+2*(n+alpha)*(n+beta)*L(k0,:))./(tmp*(1-z.^2));
        dz = L(kp1,:)./pp;
        z  = z-dz; % Newton�s formula.
        unfinished = find((abs(dz) > releps*abs(z)));
        
        if not(any(unfinished))
          break
        end
      end
      if any(unfinished)
        warning('WAFO:JRULE','too many iterations in jrule');
      end
      bp=z; %Store the root and the weight.
      wf=exp(gammaln(alpha+n)+gammaln(beta+n)-gammaln(n+1)-...
       gammaln(alpha+beta+n+1) ).*tmp*2^alfbet./(pp.*L(k0,:));
      
  end