function test_suite=test_polymul()
  initTestSuite;
end
function test_polymul_()
    x = ones (3,1); 
   y = ones (1,3); 
   b = 2; 
   c = 3; 
   assert (polymul (x, x), [1; 2; 3; 2; 1]); 
   assert (polymul (y, y), [1, 2, 3, 2, 1]); 
   assert (polymul (x, y), [1, 2, 3, 2, 1]); 
   assert (polymul (y, x), [1; 2; 3; 2; 1]); 
   assert (polymul (c, x), [3; 3; 3]); 
   assert (polymul (c, y), [3, 3, 3]); 
   assert (polymul (x, c), [3; 3; 3]); 
   assert (polymul (y, c), [3, 3, 3]); 
   assert (polymul (b, c), 6);
end
