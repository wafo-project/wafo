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
% [area,epsi] = simpson(x,f);
% [area2,epsi2] = simpson(x,area);
%
%See also: pdfnorm

%
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU Lesser General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.


%error(nargchk(1,3,nargin))
narginchk(1,3)

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
