function test_suite=test_stirlerr()
  initTestSuite;
end
function test_stirlerr_()
   assert(stirlerr(1), 0.0810614667953273, 1e-12) 
  assert(stirlerr(20), 0.00416631969199693, 1e-13) 
  assert(stirlerr(60), 0.00138887602982701, 1e-15) 
  assert(stirlerr(100), 8.33330555634921e-004, 1e-16) 
  assert(stirlerr(1000), 8.33333305555556e-005, 1e-17) 
  assert(stirlerr(realmax), 0, 1e-17)
end
