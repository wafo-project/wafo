function test_suite=test_spreading()
  initTestSuite;
end
function test_spreading_()
          % default values, respectively:  
  
   data = [10, nan, .43];  
   D = spreading(51,'cos2s',0,data); 
        % Frequency dependent direction 
   th0 = linspace(0,pi/2,257)'; 
   D = spreading(51,'cos2s',th0,data);
end
