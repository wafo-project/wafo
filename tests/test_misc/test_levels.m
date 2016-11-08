function test_suite=test_levels()
  initTestSuite;
end
function test_levels_()
  
 param = [1, 2, 3];  
 assert(levels(param), [1, 1.5, 2])
end
