function test_suite=test_rotspec()
  initTestSuite;
end
function test_rotspec_()
   S=demospec('dir'); 
  plotspec(S), hold on   
  plotspec(rotspec(S,pi/2),'r'), hold off
end
