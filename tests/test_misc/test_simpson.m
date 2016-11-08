function test_suite=test_simpson()
  initTestSuite;
end
function test_simpson_()
    x = linspace(0,4,201); 
   f = exp(-x.^2); 
   [area,epsi] = simpson(x,f); 
   assert(area, 0.886226911789523, 1e-10)
end
