%LNKRAY Link for x,F and parameters of Rayleigh distribution
%  
%   CALL  phati = lnkray(x,logR,phat,i)
%  
%     phati = parameter i as function of x, logR and phat(j) where j ~= i
%     x     = quantile
%     logR  = logarithm of the survival probability
%     
%   LNKRAY is a function connecting the quantile (x) and the survival 
%   probability (R) with the fixed distribution parameter, i.e.: 
%     phat(i) = link(x,logR,phat,i), 
%    where logR = log(Prob(X>x;phat)).
%  
%   Example % See proflog
%   
%   See also proflog



function  phati  = lnkray(x,logR,phat,ix)

phati =  x./ sqrt(-2*logR);
