function test_suite=test_strjoin()
  initTestSuite;
end
function test_strjoin_()
  
  assert(wstrjoin({'hello'}, '-'), 'hello') 
  assert(wstrjoin({'hello', 'world'}), 'hello world') 
  assert(wstrjoin({'2', '3', '4'}, '-by-'), '2-by-3-by-4') 
  assert(wstrjoin({'1', '2', '3', '4', '5'}, {' ', ',', '-', ';'}), '1 2,3-4;5') 
  assert(wstrjoin({'1', '2'}, '\n'), '1\n2') 
  assert(wstrjoin({'1'; '2'}, {'\n'}), '1\n2') 
  assert(wstrjoin({}, ','), '')
end
