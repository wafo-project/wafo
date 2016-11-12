function test_suite=test_hline()
  initTestSuite;
end
function test_hline_()
  h = hline(42,'g','The Answer') 
 
 %returns a handle to a green horizontal line on the current axes at y=42, and  
 % creates a text object on the current axes, close to the line, which reads  
 % "The Answer". 
  
 %hline also supports vector inputs to draw multiple lines at once: 
 
 hline([4 8 12],{'g','r','b'},{'l1','lab2','LABELC'}) 
 
 %draws three lines with the appropriate labels and colors. 
 close all;
end
