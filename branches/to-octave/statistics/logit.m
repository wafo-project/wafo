function z = logit(p)
%LOGIT    Logit function.
%
%	  z = lodds(p)
%
%	  The logit function is is log(p/(1-p)) and is an important part of
%	  logistic regression

%       Anders Holtsberg, 14-12-94
%       Copyright (c) Anders Holtsberg

if any(any(abs(2*p-1)>=1))
   error('A probability input please')
end
z = log(p./(1-p));
