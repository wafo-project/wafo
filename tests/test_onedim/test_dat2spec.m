function test_suite=test_dat2spec()
  initTestSuite;
end
function test_dat2spec_()
   x = load('sea.dat'); 
  S = dat2spec(x); 
  plotspec(S); 
 
  close all;
end
