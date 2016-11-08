function test_suite=test_dat2tc()
  initTestSuite;
end
function test_dat2tc_()
    x = load('sea.dat'); x1 = x(1:200,:); 
   tc = dat2tc(x1,0,'dw'); 
   plot(x1(:,1),x1(:,2),tc(:,1),tc(:,2),'ro',x1(:,1),zeros(1,200),':'); 
 
   assert(length(tc), 20); 
   close all;
end
