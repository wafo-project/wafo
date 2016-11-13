function test_suite=test_hrule()
  initTestSuite;
end
function test_hrule_()
   [x,w]= hrule(10); 
  assert(sum(x.*w), -5.2516e-019, eps); 
  assert(x, [3.43616, 2.53273, 1.75668, 1.03661, 0.34290, -0.34290,... 
             -1.03661  -1.75668  -2.53273  -3.43616], 1e-4) 
  assert(w, [7.6404e-006, 1.3436e-003, 3.3874e-002, 2.4014e-001, ... 
             6.1086e-001  6.1086e-001  2.4014e-001, 3.3874e-002,... 
             1.3436e-003  7.6404e-006], 1e-4);
end
