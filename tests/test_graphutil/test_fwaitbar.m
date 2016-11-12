function test_suite=test_fwaitbar()
  initTestSuite;
end
function test_fwaitbar_()
   h = fwaitbar(0,[],'this may take a while'); 
  for i=1:10, 
        % computation here % 
     if i==7, 
        fwaitbar(i/10,h,' soon finished') 
     else 
        fwaitbar(i/10,h) 
     end 
     pause(1) 
  end 
  close(h)
end
