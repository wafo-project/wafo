%% CHAPTERS runs selected chapter scripts and records elapsed time
function chtime = chapters(ch)
%
% CALL: chtime = chapters(ch)
% 
%   ch      = array of chapter number   
%   chtime  = array with elapsed times for the chosen chapters

chtime=zeros(1,6);
for k=ch,
    chap=['chapter' num2str(k)]
    starttid = clock;
    eval(chap)
    chtime(k) = etime(clock,starttid);
end
return

