function test_suite=test_findquot()
  initTestSuite;
end
function test_findquot_()
   t = freadtxt('findquot.m'); 
  [mq, mc] = findquot(t);  
  assert(t(mq)(2:11), 'findquot.m') 
  assert(t(mc)(1:41), '%FINDQUOT  Find quote and comment mask as')
end
