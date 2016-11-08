function test_suite=test_extralc()
  initTestSuite;
end
function test_extralc_()
    S = jonswap; 
   x = spec2sdat(S,100000,0.1,[],'random'); 
   lc = dat2lc(x); s = std(x(:,2)); 
   [lcEst,Est] = extralc(lc,s*[-2 2]); 
   [lcEst,Est] = extralc(lc,s*[-2 2],'exp,ml'); 
   [lcEst,Est] = extralc(lc,s*[-2 2],'ray,ml');
end
