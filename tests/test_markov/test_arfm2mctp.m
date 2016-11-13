function test_suite=test_arfm2mctp()
  initTestSuite;
end
function test_arfm2mctp_()
    param = [-1 1 32]; u = levels(param); 
   F = mktestmat(param,[-0.2 0.2],0.15,2); 
   F = F/sum(sum(F)); 
   Farfc = mctp2arfm({F []}); 
   F1 = arfm2mctp(Farfc); 
   cmatplot(u,u,{F+F' F1},3); 
   assert(F1(20,21:25), [0.00209800691364310, 0.00266223402503216,... 
                         0.00300934711658560, 0.00303029619424592,... 
                         0.00271822008031848], 1e-10); 
   assert(sum(sum(abs(F1-(F+F')))), 0, 1e-10) % should be zero 
 
   close all;
end
