function test_suite=test_strjoin()
  initTestSuite;
end
function test_strjoin_()
  
  assert(strjoin({'hello'}, '-'), 'hello') 
  assert(strjoin({'hello', 'world'}), 'hello world') 
  assert(strjoin({'2', '3', '4'}, '-by-'), '2-by-3-by-4') 
  assert(strjoin({'1', '2', '3', '4', '5'}, {' ', ',', '-', ';'}), '1 2,3-4;5') 
  assert(strjoin({'1', '2'}, '\n'), '1\n2') 
  assert(strjoin({'1'; '2'}, {'\n'}), '1\n2') 
  assert(strjoin({}, ','), '')
end
