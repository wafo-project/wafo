function test_suite=test_findnl()
  initTestSuite;
end
function test_findnl_()
   t = freadtxt('findnl.m'); 
  [inl,linenum] = findnl(t); 
  assert(dewhite(t(1:inl(1))), 'function  [inl,linenum] = findnl(str)') 
  assert(dewhite(t(linenum==2)), '%FINDNL Find new line characters and line numbers')
end
