function test_suite=test_dat2tp()
  initTestSuite;
end
function test_dat2tp_()
    x  = load('sea.dat'); x1 = x(1:200,:); 
   [tp, ind] = dat2tp(x1,0,'Mw'); tph = dat2tp(x1,0.3,'Mw'); 
   plot(x1(:,1),x1(:,2),tp(:,1),tp(:,2),'ro',tph(:,1),tph(:,2),'k*'); 
 
   assert(length(tp), 33); 
   assert(ind(1:5)', [12, 22, 23, 25, 27]); 
   close all;
end
