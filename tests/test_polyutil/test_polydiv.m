function test_suite=test_polydiv()
  initTestSuite;
end
function test_polydiv_()
    [q, r] = polydiv([3, 6, 9, 9], [1, 2, 3]); 
   assert(q, [3, 0]); 
   assert(r, [0, 0, 0, 9]);
end
