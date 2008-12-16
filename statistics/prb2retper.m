function retper = prb2retper(prb,npu)
%PRB2RETPER Return period from Probability of exceedance.  
%
% CALL retper = retper2prb(prb,npu)
%
%   retper  = the return period
%   prb     = probability of not exceeding level
%   npu     = the mean Number of events Per Unit time (eg Year)
%
% Example
%
% See also retper2prb

  if (any(npu) <=0)
    error('NPU must be positive!')
  end
  if (any(prb<0 | 1 <= prb))
    error('Probability must be the range [0, 1)!')
  end
  
  
  retper = 1 ./ (npu * (1 - prb));
  
