function recfig1
% RECFIG1  Location of Gullfaks C and Statfjord A platforms in The North Sea 
%

global  RECFIGNUM

if isempty(RECFIGNUM)
  disp('You must start recdemo in order to run this script')
  return
end

global map

plot(map(:,1),map(:,2))
axis([-2 12 58 64])
h = gca;
set(h,'xticklabel', [2 0:2:12]);

grid on

hold on
text( 1,62     ,'Statfjord A') 
line([1.8, 1.8], [62  ,61.2 ])
plot(1.8,61.2,'x') 
 
text(1,59.5,'Gullfaks C') 
line([1.8, 2.3 ], [59.50 ,61.20 ])
plot(2.30,61.20,'x') 

text(10.40,60.10,'Oslo') 
plot(10.80,59.85,'h')
text(10.00,63.05,'Trondheim') 
plot(10.80,63.40,'h')
text(4.00,58.80,'Stavanger') 
plot(5.52,58.90,'h')
text(3.50,60.30,'Bergen') 
plot(5.20,60.30,'h') 
hold off
wafostamp('Figure 1','(ER)')    
