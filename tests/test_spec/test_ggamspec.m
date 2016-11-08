function test_suite=test_ggamspec()
  initTestSuite;
end
function test_ggamspec_()
  wn = linspace(0,4); 
 N = 6; M = 2; 
 S = ggamspec(wn,N,M); 
 plot(wn,S) 
 assert(S(20:23),[0.705228792417578, 0.851577790974086,... 
                  0.974063777767012, 1.066965081263566], 1e-10) 
 
  close all
end
