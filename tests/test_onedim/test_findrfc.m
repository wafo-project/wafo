function test_suite=test_findrfc()
  initTestSuite;
end
function test_findrfc_()
    x = load('sea.dat');  
   tp = dat2tp(x);  
   ind = findrfc(tp(:,2),0.3);  
   waveplot(x,tp(ind,:),1,1); 
 
   assert(length(ind), 1008); 
   assert(ind(1:5)', [1, 2, 7, 8, 9]); 
   close all;
end
