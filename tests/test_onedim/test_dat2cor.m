function test_suite=test_dat2cor()
  initTestSuite;
end
function test_dat2cor_()
    x = load('sea.dat'); 
   rf = dat2cor(x,150,2); 
 
   assert(rf.R(1:3)', [ 1.0, 0.931593322326778, 0.764943054479577], 1e-10); 
   close all;
end
