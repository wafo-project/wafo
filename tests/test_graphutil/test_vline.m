function test_suite=test_vline()
  initTestSuite;
end
function test_vline_()
  h = vline(42,'g','The Answer') 
 
 % returns a handle to a green vertical line on the current axes at x=42, and  
 % creates a text object on the current axes, close to the line, which reads  
 % "The Answer". 
 
 % vline also supports vector inputs to draw multiple lines at once: 
 
 vline([4 8 12],{'g','r','b'},{'l1','lab2','LABELC'}) 
 
 % draws three lines with the appropriate labels and colors. 
 close all;
end
