function test_suite=test_hex2bin()
  initTestSuite;
end
function test_hex2bin_()
   assert(hex2bin('d8'), '11011000')
end
