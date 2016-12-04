function test_suite=test_munique()
  initTestSuite;
end
function test_munique_()
   words = strvcat('octave', 'matlab', 'scilab', 'abc', 'octave', 'matlab'); 
  [D, NR, C] = munique(words); 
  assert(D, strvcat('octave', 'matlab', 'scilab', 'abc')) 
  assert(C, [2;2;1;1]) 
  assert(D(NR,:),  words)  % i.e., D(NR,:) is exactly equal to words
end
