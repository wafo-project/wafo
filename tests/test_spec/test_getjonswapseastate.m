function test_suite=test_getjonswapseastate()
  initTestSuite;
end
function test_getjonswapseastate_()
   fetch = 10000; 
  u10   = 10; 
  ss = getjonswapseastate(u10,fetch); 
  S  = jonswap([],ss); 
  plotspec(S) 
 
  assert(ss,[0.5109333141955266, 2.7734701650790741, 2.4826805231495999,... 
            0.0753143261199741, 0.0919102507785122, 0.0162592549478214], 1e-10) 
  assert(spec2char(S,{'hm0','Tp'}), ... 
        [0.0761275095935633, 2.7588831847422388], 1e-10) 
 
  close()
end
