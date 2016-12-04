function test_suite=test_dm2unix()
  initTestSuite;
end
function test_dm2unix_()
    LF = char(10); CR = char(13); CRLF = [ CR, LF]; 
   str = cat(2,'test', CRLF, 'test2',CR,'test3',LF); 
   assert(dm2unix(str), cat(2,'test', LF, 'test2',LF,'test3',LF));
end
