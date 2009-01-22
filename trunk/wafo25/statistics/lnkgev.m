%LNKGEV Link for x,F and parameters of Generalized Extreme value distribution
%  
%   CALL  phati = lnkgev(x,logR,phat,i)
%  
%     phati = parameter i as function of x, logR and phat(j) where j ~= i
%     x     = quantile
%     logR  = logarithm of the survival probability
%     
%   LNKGEV is a function connecting the quantile (x) and the survival 
%   probability (R) with the fixed distribution parameter, i.e.: 
%     phat(i) = link(x,logR,phat,i), 
%    where logR = log(Prob(X>x;phat)).
%  
%   Example % See proflog
%   
%   See also proflog



function  phati  = lnkgev(x,logR,phat,ix)

if numel(phat)<3
  u = 0;
else
  u = phat(3);
end
 
logF = log(-expm1(logR));
sml = logR<-1;
logF(sml) = log1p(-exp(logR(sml)));
switch ix
  case 1
    error('lnkgev(x,logR,phat,i) where i=1 is not implemented!')
  case 2
    % % Reorganizing w.r.t. phat(2) (scale),
    if phat(1)~=0
      phati =  -(x-u).*phat(1)./expm1(phat(1).*log(-logF));
    else
      phati =  -(x-u)./log(-logF);
    end
  case 3
     if phat(1)~=0
      phati =  x + phat(2).*expm1(phat(1).*log(-logF))./phat(1);
    else
      phati = x+phat(2).*log(-logF);
    end
  otherwise
    error('Index to the fixed parameter is out of bounds')
end