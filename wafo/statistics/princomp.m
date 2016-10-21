%PRINCOMP Compute principal components of X
%
% CALL  [pc,z,w,Tsq] = princomp(X)
%
%   pc  the principal components
%   z   the transformed data
%   w   the eigenvalues of the covariance matrix
%   Tsq Hotelling's T^2 statistic for the transformed data
%
% Example
% x=[1,2,3;2,1,3]';
% [pc,z,w,Tsq]=princomp(x);
% m=[sqrt(2),sqrt(2);sqrt(2),-sqrt(2);-2*sqrt(2),0]/2;
% m(:,1) = m(:,1)*sign(pc(1,1));
% m(:,2) = m(:,2)*sign(pc(1,2));
%
% See also svd

% Author: Paul Kienzle
% This program is public domain.

function [pc,z,w,Tsq] = princomp(X)
  C = cov(X);
  [U,D,pc] = svd(C,0);
  if nargout>1, z = center(X)*pc; end
  if nargout>2, w = diag(D); end
  if nargout>3, Tsq = sum(stdize(z).^2,2); 
  	warning('WAFO:PRINCOMP','XXX FIXME XXX Tsq return from princomp fails some tests'); 
  end

%!shared pc,z,w,Tsq,m,x

%!test
%! x=[1,2,3;2,1,3]';
%! [pc,z,w,Tsq]=princomp(x);
%! m=[sqrt(2),sqrt(2);sqrt(2),-sqrt(2);-2*sqrt(2),0]/2;
%! m(:,1) = m(:,1)*sign(pc(1,1));
%! m(:,2) = m(:,2)*sign(pc(1,2));

%!assert(pc,m(1:2,:),10*eps);
%!assert(z,-m,10*eps);
%!assert(w,[1.5;.5],10*eps);
%!assert(Tsq,[4;4;4]/3,10*eps);

%!test
%! x=x';
%! [pc,z,w,Tsq]=princomp(x);
%! m=[sqrt(2),sqrt(2),0;-sqrt(2),sqrt(2),0;0,0,2]/2;
%! m(:,1) = m(:,1)*sign(pc(1,1));
%! m(:,2) = m(:,2)*sign(pc(1,2));
%! m(:,3) = m(:,3)*sign(pc(3,3));

%!assert(pc,m,10*eps);
%!assert(z(:,1),-m(1:2,1),10*eps);
%!assert(z(:,2:3),zeros(2),10*eps);
%!assert(w,[1;0;0],10*eps);
%!assert(Tsq,1,10*eps);
