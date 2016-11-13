function test_suite=test_wfindpeaks()
  initTestSuite;
end
function test_wfindpeaks_()
     x = (0:0.01:10);  
    S = x.^2+10*sin(3*x)+0.5*sin(50*x);  
    clf; plot(x,S); 
 % Find highest 8 peaks that are not less that 0.3*"global max" and have  
 % rainflow amplitude larger than 5. 
    assert(wfindpeaks(S',8,5,0.3), [909, 695, 482]'); 
 
   close all;
end
