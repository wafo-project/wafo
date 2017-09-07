function pdf = pdfnorm2d(X,m,S)
%PDFNORM2D Bivariate Gaussian distribution  
%
%CALL: pdf = pdfnorm2d(X,m,S)
%
%   X = 1x2 or n x 2 array of arguments
%   m = 1x2 or nn x 2 array of mean values (default zero vector)
%   S = 2x2 covariance matrix (default identity matrix)
%     If  nn=1  then  m  is filled to size  n x 2
%     else if  n ~= nn  then  X  is filled with first row to size  nn x 2
%
% Example:
% x = linspace(-5,5);
% [X1 X2] = meshgrid(x);
% f = reshape(pdfnorm2d([X1(:),X2(:)]),100,100);
% [area,epsi] = simpson(x,f);
% [area2,epsi2] = simpson(x,area);
%
% See also: pdfnorm, pdfnormnd

% History:
% Updated to Matlab 9.1 May 2017 by GL
% Modified Wafo routine pdfnorm2 by GL to accept array input  m

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
if d~=2
    error('Input X must be n x 2')
end

if nargin<2||isempty(m), m = zeros(n,d); end
if nargin<3||isempty(S), S = eye(d); end

[nn,dd]=size(m);
if nn==1,
    my=repmat(m,n,1);
elseif nn~=n,
    X=repmat(X(1,:),nn,1);
    my=m;
else
    my=m;
end

Xn = X-my;

den = (2*pi)*det(S)^(1/2);
if den< eps,
  error('Covariance matrix singular')
end

pdf = zeros(n,1);
S1 = inv(S);

for ix=1:nn
  kvf=Xn(ix,:)*S1*Xn(ix,:)';
  pdf(ix) = exp(-0.5*kvf)/den;
end
