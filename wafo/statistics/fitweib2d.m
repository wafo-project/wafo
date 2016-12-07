function phat=fitweib2d(data1,data2,varargin)
%FITWEIB2D Parameter estimates for 2D Weibull data.
% 
% CALL:  [phat,cov,pci] = fitweib2d(data1,data2,method)
%
%   phat  = maximum likelihood estimates of the parameters of the distribution
%   cov   = estimated covariance of phat
%   pci   = 95% confidense intervals for phat
%   data  = vectors of data points
%  method = 'SML' : Simultanous Maximum Likelihood estimate (Default)
%           'MML' : Marginal Maximum Likelihood estimate    
% Example:
%  sz = [10000,2];
%  R = rndweib(1,2,0,sz);
%  phat = fitweib2d(R(:,1),R(:,2),'method','sml');
%  x = linspace(0,6,200); 
%  f = pdfweib2d(x,x,phat,'mesh',true,'wdata',true);
%  mesh(f); figure(gcf+1);
%  contour(f); hold on;
%  plot(R(:,1),R(:,2),'.'); hold off;
%
%   close all;
%
% See also  fitweib, likweib2d, invnorm, hypgf, corrcoef

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


% tested on: matlab 5.2
%History:
% revised pab nov 2004
% fmins replaced with fminsearch  
% revised pab 02.11.2000
% by Per A. Brodtkorb 14.11.98

%   alpha = confidence level (default 0.05 corresponding to 95% CI)
%   g     = indices to fixed parameters not estimated     (Default [])
%   gphat = values for the fixed parameters not estimated (Default [])

% Example:
%  R = rndweib(1,2,10000,2);
%  phat0 = [2 2 0];  % set the B1 B2 and C12 respectively a priory
%  given = [2 4 5];  % = indices to the parameters phat0 set a priory
%  phat = fitweib2d(R(:,1),R(:,2),'smle',given,phat0); % estimate A1 and A2

global WAFO_WSTATS_DEFAULT_PLOTFLAG
error(nargchk(1,inf,nargin))
% Add these options?: 'shape',nan,'scale',nan,'location',0, 
options = struct('method','SML','fixpar',[],'ksign',0,...
  'alpha',0.05,'plotflag', WAFO_WSTATS_DEFAULT_PLOTFLAG,...
  'monitor',false,'optimset',optimset('disp','off','TolX',1e-6,'TolFun',1e-6,'MaxIter',600)); % default options
if (nargin==1 && nargout <= 1 && isequal(data1,'defaults'))
  phat = options; 
  return
end
options        = parseoptions(options,varargin{:});
options.method = upper(options.method);
method         = options.method;
% if (nargin < 3 ||isempty(method)), method = 'SMLE'; end
% if (nargin < 5 ||isempty(gparam)), gparam = []; end
% if (nargin < 4 ||isempty(given)),  given  = []; end
% if (nargin < 6)||isempty(alpha),   alpha  = 0.05; end


data1 = data1(:);  data2 = data2(:);
n = length(data1); n2 = length(data2);
if n~=n2,  error('data1 and data2  must have equal size'),end

rho = corrcoef(data1,data2);
rho = rho(2,1);

options.method(1) = [];
opt1 = options;
opt2 = opt1;
if ~isempty(opt1.fixpar)  
  opt1.fixpar(4:end)=[];
  opt1.fixpar(3) = 0;
  opt2.fixpar(1:2)=[];
  opt2.fixpar(3) = 0;
end

phat1 = fitweib(data1,opt1);
phat2 = fitweib(data2,opt2);
phat0 = [phat1.params(1:2),phat2.params(1:2),  sqrt(abs(rho))*sign(rho) ];

if strncmpi(method,'m',1)
  % marginal MLE
  dosearch = false;
  phat0(5)=findk(rho,phat0);  
  if ~isempty(options.fixpar)
    phat0(isfinite(options.fixpar)) = [];
  end
else
   %simultanous MLE
  dosearch = true;
end

phat = mlest(@pdfweib2d,phat0,{data1,data2},{options,'search',dosearch});
phat.dataname = inputname(1);

% if nargout > 1,
%   p_int = [alpha/2; 1-alpha/2];
%   [logL,cov] = likweib2d(phat,data1,data2,given,gparam);
%   sa  = diag(cov).';
%   pci = invnorm(repmat(p_int,1,5-length(given)),[phat; phat],[sa;sa]);
% end
 
function k=findk(rho,sparam)
  % finds k by a simple iteration scheme 
  rtol=1e-6; %relativ tolerance
  kold=sqrt(abs(rho))*sign(rho);
  kold2=kold-0.001;
  rho2=getrho(kold2,sparam);
  for ix=1:500,
    rho1=getrho(kold,sparam);
    k=kold-(rho1-rho).*(kold-kold2)./( rho1-rho2);
    %    disp(['k=' num2str(k) ]), pause
    if abs(kold-k)<abs(rtol.*k),break
    elseif abs(k)>1, tmp=sqrt(abs(rho))*sign(rho):0.000001:kold;
      [tmp2 I]=min(abs(getrho(tmp,sparam)-rho)); %linear search
      k=tmp(I);
      break,
    end
    rho2=rho1;    kold2=kold;
    kold=k;
    if ix>=500, disp('could not find k'),end
  end
  return
  
function y=getrho(k,sparam)
  % returns the correlationcoefficient based on:
 % a1=sparam(1);
  b1=sparam(2);
 % a2=sparam(3);
  b2=sparam(4);
  % and k
  if 0, % calculating correlation correct if rayleigh
    [K,E] = ellipke(k.^2); %complete elliptic integral of first and
    %second kind
    y=  (E-.5*(1-k.^2).*K-pi/4)./(1-pi/4);	
  else % alternatively and possibly slower
    qx=1./b1;qy=1./b2;
    y=  (hypgf(-qx,-qy,1,k.^2)-1).* gamma(qx).*gamma(qy)  ./ ....
	( sqrt(2.*b1.*gamma(2.*qx) - gamma(qx).^2) .* ...
	sqrt(2.*b2.*gamma(2.*qy) - gamma(qy).^2) );
    
    
  end
  return
