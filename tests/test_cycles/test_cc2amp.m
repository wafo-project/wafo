function test_suite=test_cc2amp()
  initTestSuite;
end
function test_cc2amp_()
    x=load('sea.dat');  
   TP = dat2tp(x); 
   [mM,Mm] = tp2mm(TP); 
   amp = cc2amp(mM); 
   histgrm(amp); 
 
   assert(size(amp),[1086, 1], eps); 
 
   close all;
end
