function test_suite=test_oscspec()
  initTestSuite;
end
function test_oscspec_()
   data = [0.01 4 275]; 
  S    = oscspec(data,[],2.5);  % Peak frequency at w=2.5
end
