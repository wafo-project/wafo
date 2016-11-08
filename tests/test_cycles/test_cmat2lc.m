function test_suite=test_cmat2lc()
  initTestSuite;
end
function test_cmat2lc_()
   x = load('sea.dat');  
  TP = dat2tp(x);  
  RFC = tp2rfc(TP);  
  param = [-2 2 151];  
  F = cc2cmat(param, RFC); 
  lc = cmat2lc(param, F); 
  plot(lc(:,1), lc(:,2)); 
 
  close all;
end
