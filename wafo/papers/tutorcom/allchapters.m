%% ALLCHAPTERS runs all chapter scripts and records elapsed time
%
% SCRIPT allchapters
%

chtime=zeros(1,6);

starttid=clock;
chapter1
chtime(1)=etime(clock,starttid);

starttid=clock;
chapter2
chtime(2)=etime(clock,starttid);

starttid=clock;
chapter3
chtime(3)=etime(clock,starttid);

starttid=clock;
chapter4
chtime(4)=etime(clock,starttid);

starttid=clock;
chapter5
chtime(5)=etime(clock,starttid);

starttid=clock;
chapter6
chtime(6)=etime(clock,starttid);
