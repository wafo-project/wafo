function test_suite=test_extract()
  initTestSuite;
end
function test_extract_()
   s = 'yy = foo(x1.^2,a_1*c,flag)'; 
  [names, p0, p1] = getnames(s); 
  assert(names, strvcat('yy', 'foo', 'x1', 'a_1', 'c', 'flag')) 
  rnames = strvcat('  yy', ' foo', '  x1', ' a_1', '   c', 'flag'); 
  cnames = {'yy', 'foo', 'x1', 'a_1', 'c', 'flag'}; 
  assert(extract(s, p0, p1+1), names) % Same thing  
  assert(extract(s, p0, p1+1, [], 'right'), rnames)   
  assert(extract(s, p0, p1+1, 'cell'), cnames)
end
