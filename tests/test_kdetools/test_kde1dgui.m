function test_suite=test_kde1dgui()
  initTestSuite;
end
function test_kde1dgui_()
    data = rndray(1,100,1); 
   if ismatlab,
     kde1dgui
   end
end
