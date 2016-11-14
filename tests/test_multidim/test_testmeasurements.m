function test_suite=test_testmeasurements()
  initTestSuite;
end
function test_testmeasurements_()
    type = [1 1 1]; bfs = ones(1,3);h=inf; 
   th0  = 90; 
   pos = [0 0 0;0 40 0; 20 20 0]; 
   D = testmeasurements(pos,type,th0); 
   S = dat2dspec(D,[pos type' bfs'],h); 
   plotspec(S); 
 
   close all;
end
