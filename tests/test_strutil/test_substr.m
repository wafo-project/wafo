function test_suite=test_substr()
  initTestSuite;
end
function test_substr_()
   s = 'How much wood would a Woodchuck chuck?'; 
  assert(substr(s, 0, 1), 'H')    % Get first character 
  assert(substr(s, -1, 1), '?')   % Get last character 
  assert(substr(s,  1), 'ow much wood would a Woodchuck chuck?')     % Remove first character 
  assert(substr(s,  0, -1), 'How much wood would a Woodchuck chuck') % Remove last character 
  assert(substr(s,  1, -1), 'ow much wood would a Woodchuck chuck')  % Remove first and last character
end
