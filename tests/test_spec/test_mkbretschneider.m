function test_suite=test_mkbretschneider()
  initTestSuite;
end
function test_mkbretschneider_()
   S = mkbretschneider('Hm0',6.5,'Tp' ,10);  
  fplot(S,[0,4]) 
  options = S('options'); % get options used 
  assert(fieldnames(options), {'Hm0', 'Tp', 'N', 'M', 'chkseastate'}') 
  assert(struct2cell(options), {6.5,10,5,4,'on'}') 
 
  close()
end
