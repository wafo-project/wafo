function test_suite=test_findname()
  initTestSuite;
end
function test_findname_()
   s = 'How much wood would a Woodchuck chuck?'; 
  assert(findname(s,'wo*'), strvcat('wood', 'would')) 
  assert(findname(s,'wo*','ignore'), strvcat('wood', 'would', 'Woodchuck'))
end
