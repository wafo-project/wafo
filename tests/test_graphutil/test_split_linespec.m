function test_suite=test_split_linespec()
  initTestSuite;
end
function test_split_linespec_()
   assert(split_linespec('b'), struct('color', 'b')); 
  assert(split_linespec('b--'), struct('linestyle','--', 'color', 'b')); 
  assert(split_linespec('b.'), struct('marker', '.', 'color', 'b')); 
  assert(split_linespec('b:^'), struct('linestyle',':', 'marker', '^', 'color', 'b')); 
  assert(split_linespec('blablar'), struct('color', 'r')); 
  assert(split_linespec('b', 'linecolor'), struct('linecolor', 'b'));
end
