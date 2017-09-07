function pdf = pdfnormnd(X,varargin)
%PDFNORMND Multivariate Normal probability density function.  
%
% CALL: f = pdfnormnd(X,m,S,options)
%       f = pdfnormnd(X,phat,options)
%
%        f = density function evaluated at X, size N x 1.
%        X = matrix of evaluation points, size N x D
%        m = mean, length D.               (default zeros(1,D))
%        S = Covariance matrix, size D x D (default eye(D))
%     phat = Distribution parameter struct
%            as returned from FITNORMND.  
%  options = struct with fieldnames:
%         .logp   : if TRUE, density, p, returned as log(p).
%
% Example: % Bivariate Gaussian distribution
%  x = linspace(-5,5);
%  [X1 X2] = meshgrid(x);
%  f = reshape(pdfnormnd([X1(:),X2(:)]),100,100);
%  [area,epsi] = simpson(x,f);
%  [area2,epsi2] = simpson(x,area);
%
%  assert(area2, 1.00000, 1.0e-5)
%
%See also  pdfnorm, cdfnormnd, rndnormnd, fitnormnd

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


%History
% Revised pab may 2007
% Fixed a bug thanks to Stefano Herzel <herzel@unipg.it>
% -Replaced den = (2*pi*det(S))^(d/2); with
%  den = (2*pi)^(d/2)*det(S)^(1/2);
% Revised pab 11nov2003  
% By pab 2002  
%error(nargchk(1,3,nargin))
narginchk(1,3)
options = struct('logp',false); % default options
Np = 2;
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[m,S] = deal(params{:});

[n,d]=size(X);

if isempty(m),
  m = zeros(1,d);
else
  m = m(:).';
  if d~=length(m)
    error('Mean must have length %d',d)
  end
end
if isempty(S),
  S = eye(d);
end
if any(d~=size(S))
  error('Covariance matrix must have %d dimensions',d)
end



%den = (2*pi*det(S))^(d/2);
den = (2*pi)^(d/2)*det(S)^(1/2);

if den< eps,
  error('Covariance matrix singular')
end
Xn = X-m(ones(n,1),:);

% new and fast call
if options.logp
  pdf = -0.5*sum((Xn(:,:)/S).*Xn(:,:) ,2)-log(den);
else
  pdf = exp(-0.5*sum((Xn(:,:)/S).*Xn(:,:) ,2))/den;
end
return

% old call slow
% S1 = inv(S);
% for ix=1:n
%   pdf(ix) = exp(-0.5*Xn(ix,:)*S1*(Xn(ix,:).'))/den;
% end
