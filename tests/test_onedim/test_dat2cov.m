function test_suite=test_dat2cov()
  initTestSuite;
end
function test_dat2cov_()
    x = load('sea.dat'); 
   rf = dat2cov(x,150,2); 
   assert(rf.R(1:3)', ... 
     [ 0.223686369437041, 0.208384728063068, 0.171107334682617], 1e-10); 
 
   close all;
end
