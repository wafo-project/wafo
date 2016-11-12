function test_suite=test_cc2dam()
  initTestSuite;
end
function test_cc2dam_()
    x = load('sea.dat'); TP=dat2tp(x); RFC=tp2rfc(TP);  
   bv = 3:8; 
   D = cc2dam(RFC,bv);  
   plot(bv,D,'x-'); 
 
   close all;
end
