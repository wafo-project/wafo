function test_suite=test_cmat2extralc()
  initTestSuite;
end
function test_cmat2extralc_()
    [G, Gh] = mktestmat([-1 1 64],[-0.2 0.2], 0.15,1); 
   xD = mctpsim({G Gh}, 2000); 
   Frfc = dtp2rfm(xD, 64, 'CS'); 
   [lcEst, Est] = cmat2extralc([-1 1 64], Frfc, [-0.4 0.4]); 
   lcG = cmat2lc([-1 1 64], G/sum(sum(G))); 
   lcF = cmat2lc([-1 1 64], Frfc); 
   clf; 
   semilogx(1000*lcG(:,2), lcG(:,1),... 
            lcF(:,2), lcF(:,1),... 
            lcEst.lc(:,2), lcEst.lc(:,1)); 
 
   close all;
end
