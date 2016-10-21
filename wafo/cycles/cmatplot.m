function cmatplot(in1,in2,in3,in4,in5)
%CMATPLOT Plots a cycle matrix, e.g. a rainflow matrix.
%
% CALL:  cmatplot(F)
%        cmatplot(F,method)
%        cmatplot(ux,uy,F)
%        cmatplot(ux,uy,F,method)
% 
% F      = Cycle matrix (e.g. rainflow matrix) [nxm]
% method = 1: mesh-plot (default)
%          2: surf-plot
%          3: pcolor-plot    [axis('square')]
%          4: contour-plot   [axis('square')]
%          5: TechMath-plot  [axis('square')]
%          11: From-To, mesh-plot
%          12: From-To, surf-plot
%          13: From-To, pcolor-plot   [axis('square')]
%          14: From-To, contour-plot  [axis('square')]
%          15: From-To, TechMath-plot [axis('square')]
% ux     = x-axis (default: 1:m)
% uy     = y-axis (default: 1:n)
%
% Examples:
%   param = [-1 1 64]; u=levels(param);
%   F = mktestmat(param,[-0.2 0.2],0.25,1/2);
%   cmatplot(F,1)
%   cmatplot(u,u,F,2), colorbar
%   cmatplot(u,u,F,3), colorbar
%   cmatplot(u,u,F,4)
%
% See also  cocc, plotcc

% Tested on Matlab 6.0 
%
% 1997-10-22  PJ  Korrigerat pcolor-plot. Ritar nu rätt antal rutor.
% 1998-11-03  PJ  Method 5
% Copyright (c) 1997 by Pär Johannesson
% Toolbox: Rainflow Cycles for Switching Processes V.1.0, 2-Oct-1997
% Check input arguments

ni = nargin;
no = nargout;
error(nargchk(1,5,ni));

if ni == 1     % cmatplot(F)
  F=in1;
  method=[];
  ux=[]; uy=[];
  clevels=[];
elseif ni == 2 % cmatplot(F,method)
  F=in1;
  method = in2;
  ux=[]; uy=[];
  clevels=[];
elseif ni == 3 % cmatplot(ux,uy,F)
  ux=in1;
  uy=in2;
  F=in3;
  method=[];
  clevels=[];
elseif ni == 4 % cmatplot(ux,uy,F,method)
  ux=in1;
  uy=in2;
  F=in3;
  method=in4;
  clevels=[];
elseif ni == 5 % cmatplot(ux,uy,F,method,clevels)
  ux=in1;
  uy=in2;
  F=in3;
  method=in4;
  clevels=in5;
else
  error('Wrong number of input arguments.');
end

if isempty(method)  % Default method
  method = 1;
end


% Vrid cykelmatrisen för att plotta rätt
%F = flipud(F');


% If F is a cell-array, then plot each cell in a subplot
if iscell(F)
  [N,M] = size(F);
  for i=1:N
    for j=1:M
      subplot(N,M,(i-1)*M+j)
      cmatplot(ux,uy,F{i,j},method,clevels);
    end
  end
  
else

if isempty(ux)
  ux =(1:size(F,2));  % Antalet kolumner
end
if isempty(uy)
  uy =(1:size(F,1));  % Antalet rader
end

% Make sure ux and uy are row vectors
ux = ux(:)'; 
uy = uy(:)';

n = length(F);
      
if method == 1      % mesh
  F = flipud(F');% Vrid cykelmatrisen för att plotta rätt
  mesh(ux,fliplr(uy),F)
  xlabel('min')
  ylabel('Max')
  view(-37.5-90,30)
  v = axis; axis([min(ux) max(ux) min(uy) max(uy) v(5:6)]);
elseif method == 2  % surf
  F = flipud(F');% Vrid cykelmatrisen för att plotta rätt
  surf(ux,fliplr(uy),F)
  xlabel('min')
  ylabel('Max')
  view(-37.5-90,30)
  v = axis; axis([min(ux) max(ux) min(uy) max(uy) v(5:6)]);
elseif method == 3  % pcolor
  F = flipud(F');% Vrid cykelmatrisen för att plotta rätt
  F1 = [F zeros(length(uy),1); zeros(1,length(ux)+1)];
  F2 = F1; F2(F2==0)=NaN;
  F1 = F2;
  dx=ux(2)-ux(1); dy=uy(2)-uy(1);
  ux1 = [ux ux(length(ux))+dx] - dx/2;
  uy1 = [uy uy(length(uy))+dy] - dy/2;
  pcolor(ux1,fliplr(uy1),F1)
  xlabel('min')
  ylabel('Max')
  v = axis; axis([min(ux1) max(ux1) min(uy1) max(uy1)]);
  axis('square')
elseif method == 4  % contour
  F = flipud(F');% Vrid cykelmatrisen för att plotta rätt
  if isempty(clevels)
    Fmax=max(max(F));
    clevels=Fmax*[0.005 0.01 0.02 0.05 0.1 0.2 0.4 0.6 0.8];
  end
  contour(ux,fliplr(uy),F,clevels)
  xlabel('min')
  ylabel('Max')
  v = axis; axis([min(ux) max(ux) min(uy) max(uy)]);
  axis('square')
  % Skriv höjden på konturlinjerna. (Kopia från WAT.)
%  Cstr=num2str(clevels(1),4);
%  for i=2:length(clevels)
%    Cstr=[Cstr ',' num2str(clevels(i),4)];
%  end
%  title(['ISO-lines: ' Cstr])

if 1==2
  clevels=sort(clevels);
  n_clevels=length(clevels);
  if n_clevels>12
    disp('   Only the first 12 levels will be listed in table.')
    n_clevels=12;
  end

  textstart_x=0.65;
  textstart_y=0.45;
  delta_y=1/33;
  h=figtext(textstart_x,textstart_y,'Level curves at:','norm');
  set(h,'FontWeight','Bold')

  textstart_y=textstart_y-delta_y;

  for i=1:n_clevels
    textstart_y=textstart_y-delta_y;
    figtext(textstart_x,textstart_y,num2str(clevels(i)),'norm')
  end
end % 1==2

elseif method == 5 |  method == 15 % TechMath-typ

  if isempty(clevels)
    Fmax=max(max(F));
    clevels=Fmax*[0.001 0.005 0.01 0.05 0.1 0.5 1.0];
  end
  v=clevels;
%  axis('ij');
  sym = '...x+***';
  sz = [6 20 24 8 8 8 12 16]
  
%  plot(-1,-1,sym(1),'markersize',1),hold on
  for i = 1:length(v)
    plot(-1,-1,sym(i),'markersize',sz(i)),hold on
  end
  
  for i = 1:length(v)-1
    Ind = (F>v(i)) & (F<=v(i+1));
    [I,J] = find(Ind);
%    axis('ij');
    plot(I,J,sym(i),'markersize',sz(i)),hold on
  end
  plot([1 n],[1 n],'--'), grid
  hold off

  axis([0.5 n+0.5 0.5 n+0.5])

  %legendText = sprintf('%6g < f <= %6g\n',[v(1:nv-1); v(2:nv)])
  %legendText = sprintf('<= %g\n',v(2:end))

  legendText=num2str(v(1:end)')

  legend(legendText,-1)
  
  title('From-To plot')
  xlabel('To / Standing')
  ylabel('From / Hanging')
  
  if method == 15
    axis('ij');
  end
  
elseif method == 11  % mesh
  
  mesh(ux,uy,F)
  axis('ij');
  xlabel('To')
  ylabel('From')
  view(-37.5-90,30)
  v = axis; axis([min(ux) max(ux) min(uy) max(uy) v(5:6)]);
  
elseif method == 12  % surf
  
  surf(ux,uy,F)
  axis('ij');
  xlabel('To')
  ylabel('From')
  view(-37.5-90,30)
  v = axis; axis([min(ux) max(ux) min(uy) max(uy) v(5:6)]);
  
elseif method == 13  % From-To-Matrix - pcolor
  
  F1 = [F zeros(length(uy),1); zeros(1,length(ux)+1)];
  F2 = F1; F2(F2==0)=NaN;
  F1 = F2;
  dx=ux(2)-ux(1); dy=uy(2)-uy(1);
  ux1 = [ux ux(length(ux))+dx] - dx/2;
  uy1 = [uy uy(length(uy))+dy] - dy/2;
  axis('ij');
  pcolor(ux1,uy1,F1)
  axis('ij');
  xlabel('To')
  ylabel('From')
  v = axis; axis([min(ux1) max(ux1) min(uy1) max(uy1)]);
  axis('square')
  
elseif method == 14  % contour
  if isempty(clevels)
    Fmax=max(max(F));
    clevels=Fmax*[0.005 0.01 0.02 0.05 0.1 0.2 0.4 0.6 0.8];
  end
  contour(ux,uy,F,clevels)
  axis('ij');
  xlabel('To')
  ylabel('From')
  v = axis; axis([min(ux) max(ux) min(uy) max(uy)]);
  axis('square')
  % Skriv höjden på konturlinjerna. (Kopia från WAT.)
%  Cstr=num2str(clevels(1),4);
%  for i=2:length(clevels)
%    Cstr=[Cstr ',' num2str(clevels(i),4)];
%  end
%  title(['ISO-lines: ' Cstr])

  if 1==1 % Skriv ut höjden på nivåkurverna.
    clevels=sort(clevels);
    n_clevels=length(clevels);
    if n_clevels>12
      disp('   Only the first 12 levels will be listed in table.')
      n_clevels=12;
    end

    textstart_x=0.10;
    textstart_y=0.45;
    delta_y=1/33;
    h=figtext(textstart_x,textstart_y,'Level curves at:','norm');
    set(h,'FontWeight','Bold')

    textstart_y=textstart_y-delta_y;

    for i=1:n_clevels
      textstart_y=textstart_y-delta_y;
      figtext(textstart_x,textstart_y,num2str(clevels(i)),'norm')
    end
  end

elseif method == 15  % TechMath-typ
% See: 'method == 5'

%  if isempty(clevels)
%    Fmax=max(max(F));
%    clevels=Fmax*[0.005 0.01 0.05 0.1 0.4 0.8];
%  end
%  v=clevels;
%  axis('ij');
%  sym = '...***';
%  sz = [8 12 16 8 12 16]
%  for i = 1:length(v)-1
%    Ind = (F>v(i)) & (F<=v(i+1));
%    [I,J] = find(Ind);
%    axis('ij');
%    plot(J,I,sym(i),'markersize',sz(i)),hold on
%  end
%  hold off
%
%  axis([0.5 n+0.5 0.5 n+0.5])
%
%  %legendText = sprintf('%6g < f <= %6g\n',[v(1:nv-1); v(2:nv)])
%  legendText = sprintf('<= %g\n',v(2:nv))
%
%  %legendText=num2str(v(2:nv)')
%
%  legend(legendText,-1)

end

end % if iscell(F)
