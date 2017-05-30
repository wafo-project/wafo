function test_suite=test_lexsort()
  initTestSuite;
end
function test_lexsort_()
   words = strvcat('octave', 'matlab', 'scilab', 'abc', 'Octave', 'Matlab'); 
  assert(lexsort(words), strvcat('abc','Matlab','matlab','Octave','octave','scilab'))
end
