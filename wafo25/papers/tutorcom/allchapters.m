function allchapters(sav)
%% ALLCHAPTERS runs the examples in WAFO Tutorial for Version 2.5
%
% CALL:	allchapters(sav)
%
% If  sav==1  it saves all results in  AllChapters.mat
% Default  saves=0
%

saving=0;
if nargin>0, saving=1; end

time=zeros(1,5);
starttid=clock;
chapter1
time(1)=etime(clock,starttid)
if saving, 
	save AllChapters
end 

starttid=clock;
chapter2
time(2)=etime(clock,starttid)
if saving,
	save AllChapters
end

starttid=clock;
Chapter3
time(3)=etime(clock,starttid)
if saving,
	save AllChapters
end

starttid=clock;
Chapter4
time(4)=etime(clock,starttid)
if saving,
	save AllChapters
end

starttid=clock;
Chapter5
time(5)=etime(clock,starttid)
if saving,
	save AllChapters
end
