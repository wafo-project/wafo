function phati  = lnkweib(x,logR,phat,i)
%LNKWEIB Link for x,F and parameters of the Weibull distribution
%
% CALL  phati = lnkweib(x,logR,phat,i)
%
%   phati = fixed parameter as function of x, logR and phat(j) where j ~= i
%   x     = quantile
%   logR  = logarithm of the survival probability
%   
% LNKWEIB is a function connecting the quantile (x) and the survival 
% probability (R) with the fixed distribution parameter, i.e.: 
%   phat(i) = link(x,logR,phat,i), 
%  where logR = log(Prob(X>x;phat)).
%
% Example % See proflog
% 
% See also proflog

switch i
  case 1,
    phati = (x-phat(3))./(-logR).^(1./phat(2));
  case 2,
    phati = log(-logR)./log((x-phat(3))./phat(1));
  case 3
    phati = x-phat(1).*(-logR).^(1./phat(2));
  otherwise
    error('Index to the fixed parameter is out of bounds')
end