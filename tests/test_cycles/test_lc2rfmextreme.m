function test_suite=test_lc2rfmextreme()
  initTestSuite;
end
function test_lc2rfmextreme_()
    u = (-5:0.2:5)'; lc = [u exp(-u.^2/2)]; 
   [Frfc,u,Nrfc] = lc2rfmextreme(lc); 
   cmatplot(u,u,{Frfc Nrfc},3); 
 
   close all;
end
