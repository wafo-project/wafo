function test_suite=test_binom()
  initTestSuite;
end
function test_binom_()
    b52  = binom(5,2);       % Should be 10. 
   bmax = binom(realmax,1); % Should be realmax 
   assert(b52, 10, eps); 
   assert(bmax, realmax, 1e-10);
end
