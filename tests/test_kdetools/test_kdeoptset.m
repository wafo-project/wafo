function test_suite=test_kdeoptset()
  initTestSuite;
end
function test_kdeoptset_()
   assert(kdeoptset('kdebin'), struct('kernel', 'epan', 'hs', [], ... 
         'hsMethod', 'hns', 'alpha', 0, 'L2', 1, 'inc', 128));
end
