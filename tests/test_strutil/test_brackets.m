function test_suite=test_brackets()
  initTestSuite;
end
function test_brackets_()
   s1 = '(How) much (wood (would a)) \(Woodchuck chuck\)?'; 
  [no,nc,ins] = brackets(s1,[],[],'\'); 
  
  assert(s1(find(ins)), '(How)(wood (would a))')  % find level 1 and higher enclosures 
  assert(s1(find(ins==2)), '(would a)') % find level 2 enclosures only
end
