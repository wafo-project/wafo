function test_suite=test_hex2int32()
  initTestSuite;
end
function test_hex2int32_()
   assert(hex2int32('1ab'), int32(427))
end
