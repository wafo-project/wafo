function test_suite=test_mkochihubble()
  initTestSuite;
end
function test_mkochihubble_()
   S = mkochihubble('def',2);  
  fplot(S,[0, 5]) 
 
  close()
end
