function test_suite=test_ccquad()
  initTestSuite;
end
function test_ccquad_()
  
   a=0; b=1; 
   [val,tol] = ccquad('exp(x)',a,b);  
   assert(val, exp(1)-1, tol);
end
