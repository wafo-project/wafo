function test_suite=test_perl()
  initTestSuite;
end
function test_perl_()
  
     perl -v                           % print version information 
 
     perl -le "print for @INC"         % print module search path 
     perl '-le' '"print for @INC"'     % equivalent 
     perl('-le', '"print for @INC"')   % ditto 
 
     perl -le '"print for split /;/, $ENV{PATH}"'      % print DOS path
end
