function [phat, pci,pinit]=ochi98fit(data,alpha)
%OCHI98FIT Parameter estimates and confidence intervals for Ochi data. 
%
% CALL: [phat pci] = ochi98fit(data,alpha)
%
%   phat  = [a b] = maximum likelihood estimates of the parameters of the distribution
%   pci   = 100(1-alpha) percent confidense intervals
%   data  = data matrix
%   alpha = confidence level (default 0.05 corresponding to 95% CI)
%
%
% See also  ochi98pdf 

%  Reference:
%       [1]  Michel K. Ochi,
%       "Probability distributions of peaks and troughs of non-gaussian processes"
%        Probabilistic Engineering Mechanics Vol 13 No 4 (1998) 
%       pp  291-298

% tested on: 
% history:
% revised pab nov 2004
% - replaced call to fmins with fminsearch  
% revised pab 04.11.2000
% - removed ochi98like with a call to loglike instead
% revised pab 29.02.2000
%  changed name to ochi98fit
%  Per A. Brodtkorb 14.02.99

if (nargin < 2)||isempty(alpha)
    alpha = 0.05;
end
p_int = [alpha/2; 1-alpha/2];

data1=data(:)

a = fitray(data1)*sqrt(2);
pinit=[a a];

%simultanous MLE
phat = fminsearch('loglike',pinit,[],data1,'ochi98pdf');
if nargout == 2
   [LL,cov]=loglike(phat,data1,'ochi98pdf');
   sa = diag(cov).';
   pci = invnorm(repmat(p_int,1,2),[phat; phat],[sa;sa]);
 end
 
