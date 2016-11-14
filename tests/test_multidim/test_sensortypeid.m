function test_suite=test_sensortypeid()
  initTestSuite;
end
function test_sensortypeid_()
  assert(sensortypeid(strvcat('W','v', 'n')), [12,11,1]') 
 assert(isnan(sensortypeid('rubbish')))
end
