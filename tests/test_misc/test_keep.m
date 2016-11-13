function test_suite=test_keep()
  initTestSuite;
end
function test_keep_()
  [set1data, set2data, set3data, data,test,test2,test3] = deal(1); 
 
 keep('data', 'set1data', 'set2data', 'set3data', 'test', 'test2', 'test3'); 
 assert(who(), {'data', 'set1data', 'set2data', 'set3data', 'test', ... 
               'test2', 'test3'}');  % list your workspace variables 
 
 keep('*data'); 
 assert(who(), {'data', 'set1data', 'set2data', 'set3data'}'); 
 
 keep set1data set2data set3data; 
 assert(who(), {'set1data', 'set2data', 'set3data'}');
end
