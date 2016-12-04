function test_suite=test_insert()
  initTestSuite;
end
function test_insert_()
  s1 = 'How much wood would a Woodchuck chuck?'; 
 no = [5 15]; s0 = strvcat('test','j'); 
 assert(insert(s1,s0,no), 'How testmuch wood jwould a Woodchuck chuck?') 
 assert(insert(s1,s0(1,:),no), 'How testmuch wood testwould a Woodchuck chuck?')   
 assert(insert(s1,s0,no,no+3), 'How testh wood jld a Woodchuck chuck?')   
 ix = findstr('wood',s1); 
 assert(insert(s1,s0,ix), 'How much testwood would a Woodchuck chuck?') 
 assert(insert(s1,s0,ix,ix+4), 'How much test would a Woodchuck chuck?')
end
