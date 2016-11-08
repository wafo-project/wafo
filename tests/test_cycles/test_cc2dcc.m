function test_suite=test_cc2dcc()
  initTestSuite;
end
function test_cc2dcc_()
   x = load('sea.dat'); 
  tp = dat2tp(x); 
  rfc = tp2rfc(tp); 
  param = [-2, 2, 41]; 
  dcc = cc2dcc(param,rfc); 
  u = levels(param); 
  Frfc = dcc2cmat(dcc,param(3)); 
  cmatplot(u,u,{Frfc}, 4); 
 
  close all;
end
