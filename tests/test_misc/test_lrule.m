function test_suite=test_lrule()
  initTestSuite;
end
function test_lrule_()
   [x, w] = lrule(10,2); 
  assert(sum(x.*w), 6, 1e-10) 
  [x,w]= lrule(10,2,0); 
  assert(sum(x.*w), 6, 1e-10)
end
