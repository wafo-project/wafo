function test_suite=test_mkjonswap()
  initTestSuite;
end
function test_mkjonswap_()
      options = mkjonswap('defaults'); 
     assert(fieldnames(options), {'Hm0', 'Tp', 'gamma', 'sigmaA', 'sigmaB',... 
             'Ag', 'N', 'M', 'method', 'wnc', 'chkseastate'}') 
     assert(struct2cell(options), ... 
            {7,11,[],0.07,0.09,[],5,4,'integration', 6, 'on'}') 
     options.gamma = 1; 
     S = mkjonswap(options); 
  % or alternatively 
     S = mkjonswap('gamma', 1); 
  % Plot the spectrum by 
     fplot(S,[0,5])  
  % or alternatively 
     w = linspace(0,5); 
     plot(w, S(w)) 
  
   options2 = S('options'); % get options used 
 
   S2 = mkbretschneider(options); 
   [x,y] = fplot(S,[0,4]); 
   y2 = S2(x); 
   assert(y,y2,eps) %JONSWAP with gamma=1 equals Bretscneider! 
 
   close all
end
