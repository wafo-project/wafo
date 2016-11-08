function test_suite=test_mktorsethaugen()
  initTestSuite;
end
function test_mktorsethaugen_()
   S = mktorsethaugen('Hm0',6, 'Tp',8); 
  fplot(S,[0,5]); 
  options = mktorsethaugen('defaults'); 
  assert(fieldnames(options), {'Hm0', 'Tp', 'method', 'wnc', 'chkseastate'}') 
  assert(struct2cell(options), {7,11,'integration', 6, 'on'}') 
  options = S('options'); 
  assert(fieldnames(options), {'Hm0','Tp','method','wnc','chkseastate',... 
                               'Hwoptions','Hsoptions'}') 
  assert(fieldnames(options.Hwoptions), {'Hm0', 'Tp', 'gamma', 'sigmaA', 'sigmaB',... 
             'Ag', 'N', 'M', 'method', 'wnc', 'chkseastate'}') 
  assert(fieldnames(options.Hsoptions), {'Hm0', 'Tp', 'gamma', 'sigmaA', 'sigmaB',... 
             'Ag', 'N', 'M', 'method', 'wnc', 'chkseastate'}') 
 
  close all
end
