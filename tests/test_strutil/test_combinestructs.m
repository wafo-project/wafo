function test_suite=test_combinestructs()
  initTestSuite;
end
function test_combinestructs_()
   s1 = struct('string','hello','lengths',[5 3],'test','tull');  
  s2 = struct('string','apple','lengths',[15 30],'color','red');  
  s3 = struct('string','apple','lengths',[15 30],'test','tull'); 
  s4 = struct('string','apple','lengths',[15 30],'test','tull', 'color', 'red'); 
  s5 = struct('string','hello','lengths',[5 3],'test','tull', 'color', 'red'); 
  assert(combinestructs(s1,s2,'replace'), s3) 
  assert(combinestructs(s1,s2,'merge'), s4) 
  assert(combinestructs(s1,s2,'add'), s5)
end
