function test_suite=test_dat2wa()
  initTestSuite;
end
function test_dat2wa_()
   x = load('sea.dat'); x1 = x(1:400,:);  
  [T, ind] = dat2wa(x1,0,'c2c'); % Returns crest2crest waveperiods 
  subplot(121); waveplot(x1,'-',1,1); subplot(122); histgrm(T); 
 
  assert(length(T), 22); 
  assert(ind(1:5)', [ 12, 29, 32, 40, 57]); 
  close all;
end
