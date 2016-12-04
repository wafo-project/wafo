function test_suite=test_strwcmp()
  initTestSuite;
end
function test_strwcmp_()
   s1 = 'How much wood would a Woodchuck chuck?'; 
  assert(strwcmp(s1,'*chuck*'), true) 
  assert(strwcmp(s1,'how*'), false)   
  assert(strwcmp(s1,'how*','i'), true)
end
