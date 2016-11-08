function test_suite=test_mkdspec()
  initTestSuite;
end
function test_mkdspec_()
   S = jonswap; 
  D = spreading(linspace(-pi, pi, 51), 'cos2s'); 
  Snew = mkdspec(S,D,1); 
 
  close all
end
