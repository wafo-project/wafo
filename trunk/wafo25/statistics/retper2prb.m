function prb = retper2prb(retper,npu)
%RETPER2PRB Probability of exceedance from return period. 
%
% CALL prb = retper2prb(retper,npu)
%
%   prb     =  probability of exceeding level
%   retper  = the return period
%   npu     = the mean Number of events Per unit time (eg Year)

  if (any(npu) <=0)
    error('NPU must be positive!')
  end
  if (any(retper < 1/npu))
    error('return period incompatible with NPU!')
  end
  prb =  1./(npu .* retper);
  
