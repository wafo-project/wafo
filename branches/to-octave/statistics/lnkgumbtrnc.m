%LNKGUMBTRNC Link for x,F and parameters of truncated Gumbel distribution
%  
%   CALL  phati = lnkgumbtrnc(x,logR,phat,i)
%  
%     phati = parameter i as function of x, logR and phat(j) where j ~= i
%     x     = quantile
%     logR  = logarithm of the survival probability
%     
%   LNKGUMBTRNC is a function connecting the quantile (x) and the survival 
%   probability (R) with the fixed distribution parameter, i.e.: 
%     phat(i) = link(x,logR,phat,i), 
%    where logR = log(Prob(X>x;phat)).
%  
%   Example % See proflog
%   
%   See also proflog



function  phati  = lnkgumbtrnc(x,logR,phat,ix)
% TODO % Not implemeted yet


logF = log(-expm1(logR));
sml = logR<-1;
logF(sml) = log1p(-exp(logR(sml)));

 %expba = exp(b./a);
 %tmp=exp(-exp(b(k1)./a(k1)));
 %x =-a.*log(-log(-expm1(-expba).*exp(logF) +exp(-expba)) ) + b;

switch ix
  case 1
    phati = -(x-phat(2))./log(-logF);
  case 2
    phati = x +phat(1).*log(-logF);
  otherwise
    error('Index to the fixed parameter is out of bounds')
end
error('Not implemented yet')