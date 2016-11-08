function test_suite=test_rfcfilter()
  initTestSuite;
end
function test_rfcfilter_()
       x = load('sea.dat'); 
      y = rfcfilter(x,0,1); 
            % 2. This removes all rainflow cycles with range less than 0.5. 
      y = rfcfilter(x,0.5);
end
