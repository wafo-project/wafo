function test_suite=test_spec2mom()
  initTestSuite;
end
function test_spec2mom_()
   S = demospec('dir'); 
  [m,mtext]=spec2mom(S,2,'xyt'); 
  assert(mtext, {'m0','mx','my','mt','mxx','myy','mtt','mxy','mxt','myt'}) 
  assert(m, [2.89036380451949, 1.22955597095589e-001, 2.06430283981203e-018,... 
             1.87891958511392, 6.57375173604896e-003, 8.45586322894264e-004,... 
             1.28612129886202,  1.45879204649528e-022,  8.89911740656099e-002,... 
             1.50655221045681e-018], 1e-10)
end
