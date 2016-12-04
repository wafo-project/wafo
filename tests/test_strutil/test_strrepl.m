function test_suite=test_strrepl()
  initTestSuite;
end
function test_strrepl_()
  s1 = 'How much wood would a Woodchuck chuck?'; 
 assert(strrepl(s1,'chuck','pack'), 'How much wood would a Woodchuck pack?') 
 assert(strrepl(s1,'chuck','pack','all'), ... 
                   'How much wood would a Woodpack pack?') 
 assert(strrepl(s1,'Wood', 'food','all'), ... 
                   'How much wood would a foodchuck chuck?') 
 assert(strrepl(s1,'Wood', 'food','ignore','all'), ... 
                   'How much food would a foodchuck chuck?')
end
