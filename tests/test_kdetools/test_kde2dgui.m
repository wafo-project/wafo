function test_suite=test_kde2dgui()
  initTestSuite;
end
function test_kde2dgui_()
  
  data = rndray(1,1000,2); 
  if ismatlab,
    kde2dgui
  end
end
