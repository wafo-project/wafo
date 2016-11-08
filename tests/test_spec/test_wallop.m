function test_suite=test_wallop()
  initTestSuite;
end
function test_wallop_()
    S = wallop(1.1,[6.5 10]);  
   plotspec(S) 
   
   close()
end
