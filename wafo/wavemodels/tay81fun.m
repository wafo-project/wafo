function y=tay81fun(u,x,b,pdfstr )
%TAY81FUN Internal function to tay81pdf and tay81cdf. 
%
%     Call:  y=tay81fun(u,x,b )
%
%k=find(~isreal(u))
%k=find(~isreal(x))
%k=find(~isreal(a))
%k=find(~isreal(b))
if nargin<4||isempty(pdfstr)
  pdfstr='pdf';
end
pdf=strncmp((pdfstr),'pdf',1);
%y=zeros(size(u));
[y ierr] =besselj(real(~pdf),u.*x);
switch ierr(1),
  case 0, %computation OK
  case 1, error('Illegal arguments.')
  case 2,   disp('Overflow.  Return Inf.')
  case 3,   disp('Some loss of accuracy in argument reduction.')
  case 4,   error('Complete loss of accuracy, z or nu too large.')
  case 5,  error('No convergence.  Return NaN.')
end
[tmp ierr] =besselj(0,u./sqrt(b));
switch ierr(1),
  case 0, %computation OK
  case 1, error('Illegal arguments.')
  case 2,   disp('Overflow.  Return Inf.')
  case 3,   disp('Some loss of accuracy in argument reduction.')
  case 4,   error('Complete loss of accuracy, z or nu too large.')
  case 5,  error('No convergence.  Return NaN.')
end


if pdf, % pdf
   y=y.*u.*real(tmp.^b);
else % cdf
   y=y.*u.*x.*real(tmp.^b);
end
   return



