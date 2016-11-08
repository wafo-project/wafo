function test_suite=test_freqtype()
  initTestSuite;
end
function test_freqtype_()
   S = demospec(); 
  assert(freqtype(S), 'w') 
  S2 = ttspec(S); 
  assert(freqtype(S2), 'f')
end
