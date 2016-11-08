function test_suite=test_tpextrapolate()
  initTestSuite;
end
function test_tpextrapolate_()
    x = load('sea.dat'); 
   tp = dat2tp(x,0.5); 
   tpe = tpextrapolate(tp,1,[],1); 
   clf; plot(tp(:,1),tp(:,2),'b',tpe(:,1),tpe(:,2),'r'); 
   [tpe,Pout,I] = tpextrapolate(tp,1,[],2); 
   clf;  
   plot(tpe(:,1),tpe(:,2),'b',... 
        tpe(I.min,1),tpe(I.min,2),'g*',... 
        tpe(I.max,1),tpe(I.max,2),'g*'); 
 
  close all;
end
