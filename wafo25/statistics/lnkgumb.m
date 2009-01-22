%LNKGUMB Link for x,F and parameters of Gumbel distribution
%  
%   CALL  phati = lnkgumb(x,logR,phat,i)
%  
%     phati = parameter i as function of x, logR and phat(j) where j ~= i
%     x     = quantile
%     logR  = logarithm of the survival probability
%     
%   LNKGUMB is a function connecting the quantile (x) and the survival 
%   probability (R) with the fixed distribution parameter, i.e.: 
%     phat(i) = link(x,logR,phat,i), 
%    where logR = log(Prob(X>x;phat)).
%  
%   Example % See proflog
%   
%   See also proflog



function  phati  = lnkgumb(x,logR,phat,ix)


logF = log(-expm1(logR));
sml = logR<-1;
logF(sml) = log1p(-exp(logR(sml)));
switch ix
  case 1
    phati = -(x-phat(2))./log(-logF);
  case 2
    phati = x +phat(1).*log(-logF);
  otherwise
    error('Index to the fixed parameter is out of bounds')
end