function test_suite=test_history2file()
  initTestSuite;
end
function test_history2file_()
  history2file(fullfile(waforoot,'onedim'),[],'gaussq.m')
end
