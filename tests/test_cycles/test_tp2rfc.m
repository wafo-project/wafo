function test_suite=test_tp2rfc()
  initTestSuite;
end
function test_tp2rfc_()
    x = load('sea.dat'); tp=dat2tp(x); 
   RFC1 = tp2rfc(tp);      % Default (min-to-Max cycles in residual) 
   ccplot(RFC1);  
   RFC2 = tp2rfc(tp,'CS'); % Compare with AFNOR/ASTM standard 
   [I,J] = find(RFC1(:,1)~=RFC2(:,1) | RFC1(:,2)~=RFC2(:,2)); 
   hold on; plot(RFC1(I,1),RFC1(I,2),'b+',RFC2(I,1),RFC2(I,2),'rx'); hold off; 
 
   close all;
end
