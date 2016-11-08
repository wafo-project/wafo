function test_suite=test_mm2lc()
  initTestSuite;
end
function test_mm2lc_()
    x = load('sea.dat'); tp = dat2tp(x); mM = tp2mm(tp); 
   lc = mm2lc(mM);
end
