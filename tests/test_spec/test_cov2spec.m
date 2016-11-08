function test_suite=test_cov2spec()
  initTestSuite;
end
function test_cov2spec_()
   R = createcov; 
  L = 129;   
  R.t = linspace(0,75,L).'; 
  R.R = zeros(L,1); 
  win = parzen(40);   
  R.R(1:20) = win(21:40); 
  S0 = cov2spec(R); 
  S = jonswap; 
  S1 = cov2spec(spec2cov(S)); 
  assert(all(abs(S1.S-S.S)<1e-4) ,'COV2SPEC')
end
