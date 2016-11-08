function test_suite=test_tp2mm()
  initTestSuite;
end
function test_tp2mm_()
    x = load('sea.dat'); 
   TP = dat2tp(x); 
   [mM,Mm] = tp2mm(TP); 
   ccplot(mM);
end
