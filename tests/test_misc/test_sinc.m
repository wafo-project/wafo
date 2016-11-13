function test_suite=test_sinc()
  initTestSuite;
end
function test_sinc_()
   x =linspace(-5,5)';  
  plot(x,sinc(x)); 
 
  close();
end
