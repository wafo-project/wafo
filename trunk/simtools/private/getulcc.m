function [row,col,num]=getulcc(M)
%GETULCC Finds the upper left element > 0 in a matrix.

%  Does not work if several element lay on a diagonal line and the upper left
%  element is not unique.

%  Copyright 1993, Mats Frendahl, Dept. of Math. Stat., University of Lund.

n=length(M);
for i=n-1:-1:0
  index=find(diag(fliplr(M),i)>0);
  if ~isempty(index)
    row=index;
    col=(n-i)-(row-1);
    num=M(row,col);
    return
  end
end

