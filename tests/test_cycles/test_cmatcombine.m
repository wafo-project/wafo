function test_suite=test_cmatcombine()
  initTestSuite;
end
function test_cmatcombine_()
    F1 = triu(ones(8),0); 
   F2 = 2*F1; 
   [F, Lim, FF1, FF2] = cmatcombine(F1, F2, 2); 
 
   assert(Lim, struct('range', 2, 'min',  8,'max', 1)); 
   assert(F, FF1+FF2, eps); 
 
   Lim.range=2; Lim.min=4; Lim.max=4; 
   [F, Lim, FF1, FF2] = cmatcombine(F1,F2,Lim); 
 
   assert(Lim, struct('range', 2, 'min',  4,'max', 4)); 
   assert(F, FF1+FF2, eps);
end
