function test_suite=test_grule()
  initTestSuite;
end
function test_grule_()
   [x,w] = grule(11,0,3); 
   assert(sum(exp(x).*w), 19.0855369231877, 1e-10);
end
