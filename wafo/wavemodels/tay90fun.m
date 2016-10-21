function y=tay90fun(u,x,b)
%TAY90FUN Internal integrand function to tay90pdf 
%
%  Call:  y = tay90fun(u,x,b)
%
%    u = integration variable
%    x = normalized waveheight
%    b = correlation coefficient, i.e.,  cov(At^2,Ac^2)
%    

%k=find(~isreal(u))
%k=find(~isreal(x))
%k=find(~isreal(a))
%k=find(~isreal(b))
% if nargin<4|isempty(pdfstr)
%   pdfstr='pdf';
% end
%pdf=strcmp(pdfstr,'pdf');
A=0.25./(1-b.^2);
[y ierr] =besseli(0,u.*x.^2.*b.*A,1);
switch ierr(1),
  case 0, %computation OK
  case 1, error('Illegal arguments.')
  case 2, error('Overflow.  Return Inf.')
  case 3, disp('Some loss of accuracy in argument reduction.')
  case 4, error('Complete loss of accuracy, z or nu too large.')
  case 5, error('No convergence.  Return NaN.')
end


y=sqrt(u)/2.*x.^3.*A.*exp(-A.*x.^2.*(2-u-abs(u.*b)) ).*y;

%y=sqrt(u)/2.*x.^3./(1-b.^2).*exp(-x.^2./(1-b.^2).*(1-u/2-abs(real(u.*b/2)))).*y;%/0.73263369802690;
   return
