function test_suite=test_humps()
  initTestSuite;
end
function test_humps_()
     x = linspace(0,1); 
    y = humps(x); 
    h = plot(x,y); 
        
    close();
end
