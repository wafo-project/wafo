function test_suite=test_hex2int16()
  initTestSuite;
end
function test_hex2int16_()
   assert(hex2int16('1ab'), int16(427))
end
