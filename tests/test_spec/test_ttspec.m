function test_suite=test_ttspec()
  initTestSuite;
end
function test_ttspec_()
            % angle in radians to angle in degrees.   
    S   = demospec('dir'); 
    Sf  = ttspec(S); 
    Sf1 = ttspec(Sf,'f','d'); % = ttspec(Sf,'f','degrees');  
    plotspec(S,3), figure(2) 
    plotspec(Sf,3), figure(3) 
    plotspec(Sf1,3) 
 
    close('all')
end
