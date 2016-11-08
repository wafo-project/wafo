function test_suite=test_qrule()
  initTestSuite;
end
function test_qrule_()
    [bp,wf] = qrule(10); 
   assert(sum(bp.^2.*wf), 0.666666666666666, 1e-10) % integral of x^2 from a = -1 to b = 1 
   [bp,wf] = qrule(10,2); 
   assert(sum(bp.^2.*wf), 0.886226925452758, 1e-10)  % integral of exp(-x.^2)*x.^2 from a = -inf to b = inf 
   [bp,wf] = qrule(10,4,1,2); 
   assert(sum(bp.*wf), 0.266666666666668, 1e-10)   % integral of (x+1)*(1-x)^2 from  a = -1 to b = 1
end
