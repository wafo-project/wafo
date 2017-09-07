% ALLCHAPTERS testfile for all chapters 

tid=zeros(1,6);
starttid=clock;
Chapter1
tid(1)=etime(clock,starttid);

starttid=clock;
Chapter2
tid(2)=etime(clock,starttid);

starttid=clock;
Chapter3
tid(3)=etime(clock,starttid);

starttid=clock;
Chapter4
tid(4)=etime(clock,starttid);

starttid=clock;
Chapter5
tid(5)=etime(clock,starttid);

starttid=clock;
Chapter6
tid(6)=etime(clock,starttid);

disp(['ElapsedTime  ' num2str(tid(1)) ', ' num2str(tid(2)) ', ' num2str(tid(3))...
    ', ' num2str(tid(4)) ', ' num2str(tid(5)) ', ' num2str(tid(6)) ' sec'])