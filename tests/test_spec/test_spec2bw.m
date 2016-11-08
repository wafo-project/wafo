function test_suite=test_spec2bw()
  initTestSuite;
end
function test_spec2bw_()
    S=demospec; 
   assert(spec2bw(S,[1 2 3 4]), [0.895608469201822, 0.230162951768626,... 
                                  0.444843196973911, 3.002022823753302], 1e-10)
end
