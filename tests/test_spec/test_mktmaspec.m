function test_suite=test_mktmaspec()
  initTestSuite;
end
function test_mktmaspec_()
     options = mktmaspec('defaults'); % default options 
    options.h = 20;options.gamma = 1; 
    S = mktmaspec(options);   % Bretschneider spectrum Hm0=7, Tp=11 
    for h = [10 21 42] 
      S = mktmaspec('h',h); 
      fplot(S,[0 3]), hold on 
    end
end
