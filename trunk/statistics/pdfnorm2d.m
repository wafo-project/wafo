function pdf = pdfnorm2d(X,m,S)
%PDFNORM2D Bivariate Gaussian distribution  
%
% CALL: pdf = pdfnorm2d(X,m,S)
%
%   X = 
%   m = mean              (default zero vector)
%   S = Covariance matrix (default identity matrix)
%
% Example:
% x = linspace(-5,5);
% [X1 X2] = meshgrid(x);
% f = reshape(pdfnorm2d([X1(:),X2(:)]),100,100);
% [area,epsi] = boole(x,f);
% [area2,epsi2] = boole(x,area);
%
%See also: pdfnorm

error(nargchk(1,3,nargin))


[n,d]=size(X);

if nargin<2||isempty(m), m = zeros(1,d); end
if nargin<3||isempty(S), S = eye(d); end


Xn = X-m(ones(n,1),:);

den = (2*pi*det(S))^(d/2);
if den< eps,
  error('Covariance matrix singular')
end

pdf = zeros(n,1);
S1 = inv(S);
for ix=1:n
  pdf(ix) = exp(-0.5*Xn(ix,:)*S1*(Xn(ix,:).'))/den;
end
