function test_suite=test_cocc()
  initTestSuite;
end
function test_cocc_()
    x=load('sea.dat'); 
   tp = dat2tp(x); 
   rfc = tp2rfc(tp,'CS');           % Rainflow cycles 
   RFM = dat2rfm(x,0,[-3 3 100]);  % Rainflow matrix 
   ERFM = rfmextrapolate(RFM);     % Extrapoalted RFM 
   cocc([-3 3 100],rfc,RFM); 
   cocc([-3 3 100],rfc,ERFM);
end
