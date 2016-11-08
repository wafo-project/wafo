function test_suite=test_bretschneider()
  initTestSuite;
end
function test_bretschneider_()
   S = bretschneider(1.5,[6.5 10]); 
  assert(S.S(100:105)', [ 5.60674057893602, 5.70599250300968, ... 
                            5.79072538638665, 5.86131067761872, ... 
                            5.91820524913600, 5.96193576151422], 1e-9) 
  assert(S.w(100:105)', [0.580078125000000, 0.585937500000000, ... 
                        0.591796875000000, 0.597656250000000, ... 
                        0.603515625000000, 0.609375000000000], 1e-9) 
  plotspec(S) 
 
  close all
end
