function [M,V] = pdfstat(f,condon,cvar)
%PDFSTAT  Mean and variance for the PDF 2D distribution
% 
%   CALL:  [M,V] = pdfstat(phat,condon,cvar)
%     
%     M,V   = mean and variance, respectively
%     pdf   = pdf structure array (see createpdf)
%     condon = 0 returns marginal mean and variance for X1, X2 (default)
%              1 returns conditional mean and variance of X2 given X1 
%              2 returns conditional mean and variance of X1 given X2
%      cvar  = conditional variable, i.e.,x1 or x2 depending on condon.
% 
%  Example:
%   x1=linspace(0,6)';
%   par = {2 2  3 2.5 .8};
%   [m, v] = momweib2d(par{:}); %  mean and covariance
%   [m, v] = momweib2d(par{:},'condon',2,'cvar',x1);
%   plot(x1,m,'r--', x1,sqrt(v),'g-'); % conditional mean and standard deviation.
% 
%   [z1,z2] = rndweib2d(par{:},500,1);
%   f = kdebin([z1,z2]);
%   [M1,V1]=pdfstat(f,2,x1);
%
%   plot(x1,M1,'r--',x1,sqrt(V1),'k-');
%   title(' Conditional mean and standard deviation');
%   legend('E(x1|x2)','std(x1|x2)');
%   xlabel('x2');
%
%   close all;
%  % Does not work correctly yet!
% 
%  See also  kdebin, pdfplot
  
  
  
%tested on: matlab 5.2
% history:
%  by Per A. Brodtkorb July 2004

% TODO % needs further testing.
if nargin<3||isempty(cvar) ,
  cvar=[];
end
  
if (nargin <2) ||  isempty(condon), 
 condon=0;
end
if (condon~=0 )&&(nargin ==0), 
  error('Requires one input argument the levels to condition on.'); 
end

Ndim = length(f.x);
if (Ndim>2 || Ndim<=0)
  error('Wrong dimension')
end
f0 = 1;
switch condon, %marginal stats
 case 0,
  switch Ndim
   case {1},
    M =  trapz(f.x{1},f.x{1}(:).*f.f(:));
    V =  trapz(f.x{1},f.x{1}(:).^2.*f.f(:))-M.^2;
   case {2},
    pdf1 = trapz(f.x{2},f.f).';
    pdf2 = trapz(f.x{1},f.f.').';
    M(1) =  trapz(f.x{1},f.x{1}(:).*pdf1);
    V(1) =  trapz(f.x{1},f.x{1}(:).^2.*pdf1)-M(1).^2;
    M(2) =  trapz(f.x{2},f.x{2}(:).*pdf2);
    V(2) =  trapz(f.x{2},f.x{2}(:).^2.*pdf2)-M(2).^2;
  end
 case 1 , % conditional stats of X2 given X1
  %f0 = simpson(f.x{1},simpson(f.x{2},f.f));
  xi = f.x{1}(:);
  XI2 =  repmat(f.x{2}(:),1,length(xi));
  M = trapz(f.x{2},XI2.*f.f/f0).';
  V = trapz(f.x{2},XI2.^2.*f.f/f0).' - M.^2;
 case 2, % conditional stats of X1 given X2
  %f0 = simpson(f.x{1},simpson(f.x{2},f.f));
   
  xi=f.x{2}(:);
  XI2 = repmat(f.x{1}(:),1,length(xi));
  M = trapz(f.x{1},XI2.*f.f.'/f0).';
  V = trapz(f.x{1},XI2.^2.*f.f.'/f0).'-M.^2;
end
if (~isempty(cvar) && Ndim==2),
  
  M = interp1(xi,M,cvar);
  V = interp1(xi,V,cvar);
end

return



