function xstar=mindist(g,b,u,n0,epsi)
% MINDIST Finds minimal distance to the origin on the surface b'*x+x'*diag(g)*x=u
% 
%   CALL: xstar=mindist(g,b,u,n0,eps);
% 
%       u  = levels
%       n0 = number of non zero elements in the starting point
%       eps = accuracy in the iteration
%
% Iterative solution to the problem om finding the point of minimal
% distance to origin on the surface b'*x+x'*diag(g)*x=u, where x0
% is the starting value.


if nargin<5
   epsi=1e-12;
end
niter=1000;
iter=1;
n1=length(g);
n0=min(n1,n0);
xstar=mindist3(g(end-n0+1:end),b(end-n0+1:end),u);
if n0<n1
    xn=[zeros(n1-n0,1);xstar];
b=b(:);
g=g(:);
xnold=xn+1;
while iter<niter && sum(abs(xn-xnold))>epsi
   a=1-b'*xn/u-g'*xn.^2/u;
   A=-b/u-2*g.*xn/u;
   xnold=xn;
   xn=(A'*xn-a)/(A'*A)*A;
   iter=iter+1;
end
xstar=xn;
if iter==niter
   warning(['The algorithm didn''t converged after ' int2str(niter) 'times for level u=' int2str(u)])
end
end



