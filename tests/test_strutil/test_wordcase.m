function test_suite=test_wordcase()
  initTestSuite;
end
function test_wordcase_()
  s1 = 'How much wood would a Woodchuck chuck?'; 
 assert(wordcase(s1), 'How Much Wood Would A Woodchuck Chuck?')
end
