function test_suite=test_geth1line()
  initTestSuite;
end
function test_geth1line_()
   assert(geth1line('geth1line',1), ... 
   'geth1line     - Extracts the first comment line (the H1 line) of a m-file'); 
  assert(geth1line('geth1line',0), ... 
    'GETH1LINE     Extracts the first comment line (the H1 line) of a m-file'); 
  assert(geth1line('geth1line',1,7), ... 
    'geth1line - Extracts the first comment line (the H1 line) of a m-file');
end
