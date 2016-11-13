function test_suite=test_mctp2arfm()
  initTestSuite;
end
function test_mctp2arfm_()
    param = [-1 1 32]; u = levels(param); 
   F = mktestmat(param,[-0.2 0.2],0.15,2); 
   Frfc = mctp2rfm({F []}); 
   Farfc = mctp2arfm({F []}); 
   cmatplot(u,u,{Frfc Farfc},3); 
   assert(sum(sum(abs((Farfc+Farfc')-(Frfc+Frfc')))), 0, 1e-10) % should be zero 
 
  close all;
end
