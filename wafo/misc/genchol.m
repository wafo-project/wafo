function [L,P,r] = genchol(A,tol)
%GENCHOL Generalized Cholesky factorization
% 
% CALL: [L,P,r] = genchol(A,tol);
% 
%  L  = lower triangular matrix, Cholesky factor
%  P  = permutation vector
%  r  = matrix rank, i.e., the number of eigenvalues larger than tol.
%       This is an estimate of the number of linearly
%       independent rows or columns of a matrix A.
%  A  = symmetric and semi-positive definite matrix
% tol = tolerance used in factorization, i.e., norm(A(P,P)- L*L.') <= tol  
%
% GENCHOL computes the generalized Cholesky decomposition,
%          
%   A(P,P) = L*L.' 
% 
% where L is lower triangular and P is a permutation vector.
%
% Example
%  H = hilb(10);
%  tol   = 1e-6;
%  [L,P] = genchol(H,tol);
%  spy(L*L.'-H(P,P));
%
% close all;
% 
% See also  chol  
  
%error(nargoutchk(2,2,,nargout))

[m,n] = size(A); 
if (m ~= n),   error('input matrix must be square'); end
if nargin<2||isempty(tol),  
  tol = 1e-12;
end
%tol = tol/(100*n*n);%(4*sqrt(n));



L = tril(A); 
%D = diag(A); 
%Dmax = max(D);
%tol = Dmax*tol;
%localTolerance = max(eps*Dmax,tol*1e-3);
localTolerance = 0;%tol/(20*n^2);

P = 1:n; % permuation vector
%x = P;

k = 1;
nullity = 0;
while (k<=n)
   % Find next pivot
   D = diag(L,-nullity);
   k0 = k-nullity;
	n0 = n-nullity;
   [big,imax] = max(D(k0:n0));
   
   if (big>tol), %D(imax)>tol)
      imax = imax+k-1;
      if imax~=k
         % Swap the new pivot at imax with k
         P([imax,k]) = P([k,imax]);% [P(imax),P(k)] = deal(P(k),P(imax));
         L = rowcolChange(L,k,imax,n,nullity);     
      end
      L(k,k0) = sqrt(L(k,k0)); % STD(Xk|X1,X2,...,Xk-1)
      k1 = k;
      for i = k+1:n,
         %disp(i)
         %tmp = 0;
         %for j = 1:k-1, 
         %Cov(Xi,Xj|X1,X2,...Xk-2)*Cov(Xj,Xk|X1,X2,...Xk-2)/Cov()
         %tmp = tmp + L(i,j) * L(k,j); 
         %end
         %L(i,k) = L(i,k)-tmp; % Cov(Xi,Xk|X1,X2,...Xk-1)
         %L(i,k) = L(i,k) / L(k,k);            
         
         j =  1:k0-1;
         % Cov(Xi,Xk|X1,X2,...Xk-1)/STD(Xk|X1,X2,...Xk-1)
         L(i,k0) = (L(i,k0)-sum(L(i,j).*L(k1,j)))/L(k1,k0);
         % Var(Xi|X1,X2,...Xk)
         i0 = i-nullity;
         L(i,i0) = L(i,i0) - L(i,k0)*L(i,k0);        
         
         if (L(i,i0)<=localTolerance) % && norm(L(nnDet,k+1:nnDet),'inf')<=tol))
            %if (D(i)<=-sqrt(tol)) % make sure we are not too restrictive
            %   warning('Matrix is not positive semi-definite!')
            %end
     
            if (k+1<i)
               % Swap the singular pivot at i with k+1
               P([i,k+1]) = P([k+1,i]);
               L = rowcolChange(L,k+1,i,n,nullity);  
            end
            % shift
            if nullity>0
               n0 = n-nullity;
               L(k+1:n,k0+1:n0) = L(k+1:n,k0+2:n0+1); 
            else
               L(k+1:n,k+1:n-1) = L(k+1:n,k+2:n);
               L(n,n) = 0;
            end
            nullity = nullity+1;
            k = k+1;
         end
    
      end      
   else
      if (big<=-sqrt(tol))
         warning('Matrix is not positive semi-definite!')
      end
      k0 = k-nullity;
      L(k:end,k0:end) = 0;
      nullity = n-k0+1;
      break;
   end
   k = k+1;
end  %while k
r = n-nullity;
%flop
%nDet
return



function L = rowcolChange(L,i,j,n,nullity)
%rowcolChange exchange column and rows of i and j,
%but only the lower triangular part

if (i<j)
   k = i;
   imax = j;
else
   k = j;
   imax = i;
end	
   
   %for 
   k0 = k-nullity;
   if 1<k0
       iz=1:k0-1;
      [L(k,iz),L(imax,iz)] = deal(L(imax,iz),L(k,iz));
   end	
   
   imax0 = imax-nullity;
   [L(k,k0),L(imax,imax0)] = deal(L(imax,imax0),L(k,k0));
   %for
   iz=k+1:imax-1;
   iz0= k0+1:imax0-1;
   if any(iz)
      [L(iz,k0),L(imax,iz0)] = deal(L(imax,iz0).',L(iz,k0).');
   end
%   for 
iz=imax+1:n;
if any(iz)
  [L(iz,k0),L(iz,imax0)] = deal(L(iz,imax0),L(iz,k0));
end
