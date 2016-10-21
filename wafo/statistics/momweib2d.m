function [M ,V ,Tm ,cvar, tolm, tolv] = momweib2d(varargin)
% MOMWEIB2D  Mean and variance for the 2D Weibull distribution 
%
%  CALL:  [M,V,Tm] = momweib2d( a1,b1,a2,b2,c12,options ) 
%         [M,V,Tm] = momweib2d(phat,options ) 
%
%    M , V, Tm  = mean, variance and modal value, respectively
%   ai, bi, c12 = parameters of the distribution
%     phat = Distribution parameter struct
%            as returned from FITWEIB2D.  
%  options = struct with fieldnames:
%     .logp    : if TRUE, probability, p, returned as log(p).
%     .condon  : 0 Return mean, covariance and modal value of X1 and X2 (default)
%                1 Return the conditional values given X1 (cvar)
%                2 Return the conditional values given X2 (cvar)
%     .cvar    : vector of conditional values 
%                 (default depending on marginal mean and variance)
%     .releps  : Requested relative error (for the integrations) (Default 1e-3) 
%
% The distribution is defined by its PDF:
%    f(X1,X2)=B1*B2*xn1^(B1-1)*xn2^(B2-1)/A1/B1/N*...
%             exp{-[xn1^B1 +xn2^B2 ]/N }*I0(2*C12*xn1^(B1/2)*xn2^(B2/2)/N) 
%  where 
%    N=1-C12^2, xn1=X1/A1,  xn2=X2/A2 and 
%    I0 is the modified bessel function of zeroth order.
%
% (The marginal distribution of X1 is weibull with parameters A1 and B1) 
% (The marginal distribution of X2 is weibull with parameters A2 and B2) 
%
% Examples:
%  par = {2 2  3 2.5 .8};
%   [m v] = momweib2d(par{:}) %  mean and covariance
%   x = linspace(0,6)';
%   [m v] = momweib2d(par{:},'condon',2,'cvar',x);
%   plot(x,m,'r--', x,sqrt(v),'g-') % conditional mean and standard deviation.
%
% See also pdfweib2d, cdfweib2d, rndweib2d, fitweib2d, ellipke, hypgf, gamma

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



%   References:
%     Dag Myrhaug & Håvard Rue
%  Journal of Ship Research, Vol 42, No3, Sept 1998, pp 199-205 

%tested on: matlab 5.1
% history:
%  by Per A. Brodtkorb 13.11.98


error(nargchk(3,15,nargin))
Np = 5;
options = struct('condon',0,'cvar',[],'releps',1e-4); % default options
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[a1,b1,a2 b2,c12] = deal(params{:});

if isempty(a1)||isempty(b1)||isempty(a2)||isempty(b2)||isempty(c12)
  error('Requires either 7 input arguments or that input argument 3 is FDATA.'); 
end
[iscmn  a1 b1 a2 b2 c12] = iscomnsize(a1,b1,a2,b2,c12);
if ~iscmn
    error('Requires non-scalar arguments to match in size.');
end

[m1,v1]=momweib(a1,b1); % 
[m2,v2]=momweib(a2,b2);
tm1=zeros(size(a1));tm2=zeros(size(a2));


k=find(b1>1);
if any(k),
  tm1(k)=a1(k).*(1-1./b1(k)).^(1./b1(k));
end
k2=find(b2>1);
if any(k2),
  tm2(k2)=a2(k2).*(1-1./b2(k2)).^(1./b2(k2));
end
cvar = options.cvar;
tol   = options.releps;
switch options.condon
  case 0,    % mean and covariance
    covar=zeros(size(a2));
    ok = a1 > 0 & b1 > 0 &a2 > 0 & b2 > 0 | abs(c12)<1;
    k = find(~ok);
    if any(k)
      tmp   = NaN;
      covar(k) = tmp(ones(size(k))); 
    end
    k = find(ok);
    if any(k),
      if 0, % calculating covar correct for Rayleigh
        [K,E] = ellipke(c12(k).^2); %complete elliptic integral of first and
        %second kind
        covar(k) =  sqrt(v1(k).*v2(k)).*(E-.5*(1-c12(k).^2).*K-pi/4)./(1-pi/4);
      else % calculating covar correct for Weibull
        qx=1./b1(k);qy=1./b2(k);
        covar(k)=   sqrt(v1(k).*v2(k)).*(hypgf(-qx,-qy,1,c12(k).^2)-1).* gamma(qx).*gamma(qy)./ ....
              ( sqrt(2.*b1(k).*gamma(2*qx) - gamma(qx).^2) .* ...
                sqrt(2.*b2(k).*gamma(2*qy) - gamma(qy).^2) );	
      end
    end
    M=[m1 m2];
    V=[v1 covar;covar' v2 ];
    Tm=[tm1,tm2];
  case {1,2}, 
    sparam={a1 b1 a2 b2 c12};
   
    if options.condon==1,%conditional mean and variance given x1  
      txt='x2';      vr=v2;      mn=m2; mc=m1;
      porder=[3:4 1:2 5];
   else%conditional mean and variance given x2  
     txt='x1'; vr=v1; mn=m1; mc=m2;
     porder=1:5;
   end
     if isempty(cvar),cvar=linspace(0 ,mn+3*sqrt(vr),30)'; end
    if 1
      
       xinf=mn+mn./mc.*cvar   +15*sqrt(vr); % infinite value for x1 or x2
      %tmp=input(['Enter an infinite value for ', txt, ':  (' , num2str(xinf), ')']);
      %if ~isempty(tmp), xinf=tmp;end
      %disp(['Infinite value for ' ,txt,' is set to ' num2str(xinf(:)')])
    end
    
    M=zeros(size(cvar));Tm=M;  
    %do the integration with a quadrature formulae
     [M   tolm]=gaussq(@pdfweib2d,0,xinf,tol,[],cvar,sparam{porder},'condon',3);  %E(x2|x1) or E(x1|x2)
     [V  tolv]=gaussq(@pdfweib2d,0,xinf,tol,[],cvar,sparam{porder},'condon', 4);%E(x2^2|x1) or E(x1^2|x2) 
     V=V-M.^2;%var(x2|x1) or var(x1|x2)
     k=find(xinf<M+6*sqrt(V));
     if any(k),
       %xinf(k)=M(k)+6*sqrt(V(k))+xinf(k); % update the infinite value                                           
       
       
       disp(['Changed Infinite value for ', txt]);%, ' to ',   num2str(xinf(k))])
     end
     
     if sparam{porder(2)}>1 && nargout>2
       %modalvalue
       Nint=length(cvar(:));
       for ix=1:Nint,
         Tm(ix)=fminbnd('-pdfweib2d(x,P1,P2,P3)',max(0,M(ix)-sqrt(V(ix))),...
           M(ix)+sqrt(V(ix)),[],cvar(ix),sparam{porder},2 );
       end
     end     
 end
 
 
 % [M(ix)   tolm(ix)]=quadg('pdfweib2d',0,xinf(ix),tol,[],cvar(ix),sparam(porder),3);  %E(x2|x1) or E(x1|x2)
 % [V(ix)  tolv(ix)]=quadg('pdfweib2d',0,xinf(ix),tol,[],cvar(ix),sparam(porder),4);%E(x2^2|x1) or E(x1^2|x2) 
 %V(ix)=V(ix)-M(ix)^2; %var(x2|x1) or var(x1|x2)
 %if Nint>ix & xinf(ix)<M(ix)+6*sqrt(V(ix)),
 %  xinf(ix+1)=xinf(ix+1)+M(ix)+10*sqrt(V(ix))-xinf(ix); % update the infinite value                                           
 %  disp(['Changed Infinite value for ', txt, ' to ',   num2str(xinf(ix))])
 %end