function test_suite=test_rfmextrapolate()
  initTestSuite;
end
function test_rfmextrapolate_()
    [G,Gh] = mktestmat([-1 1 64],[-0.2 0.2], 0.15,1); 
   xD = mctpsim({G Gh},2000); 
   Frfc = dtp2rfm(xD,64,'CS'); 
   Fest = rfmextrapolate(Frfc,[],1); 
   Grfc = mctp2rfm({G Gh}); 
   cmatplot({Frfc Fest; Grfc G},4); 
 
   close all;
end
