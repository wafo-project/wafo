function Mv=seamovie(Y,s,Wavename)
%SEAMOVIE Makes a movie of simulated sea and optionally saves it
%  
%CALL: Mv = seamovie(Y,s,Wavename)
%  
%      Mv = movie
%
%      Y = struct with 2d or 3d simulation (from spec2wave or spec2field)
%      s = type of plot if 3d: 
%               if s=1 then surf-plot, (default)
%               if s=2 contour,
%               else gray-scale view 
%
%      Wavename = optional namestring for avi-file ('MyWave.avi') 
%              If absent, no avi-file is produced
%              If MyWave.avi exists in working folder 
%                   a new random name is given
%              The avi-option uses VideoWriter. 
%                   movie2avi works for Matlab ver 8.5 
%                   but not for ver 8.1
%  
% See also  spec2field, spec2wave, movie, getframe

% Tested on Matlab 8.6, 8.5, 9.1 for avi-option  
% Tested on Matlab 8.1, 5.3 - avi-option does not work correctly 
%
% History
% revised by GL 18-June-2017 to save wave movie using VideoWriter
% revised by GL March 2015 to use WafoL axis convention and plotting
% revised pab June 2005
% -fixed a bug: in matlab7: "Mv =[]; Mv(j) = getframe;" 
%       does not work, now fixed.
% revised es 20.06.00 if wrong dimension then message and return, 
%       not error
% Revised by es 13.06.00 more dimension checks  
% By es 23.05.00

if nargin>2,
    if verLessThan('matlab','8.5'),
        saveavi=0;
        disp('Too old Matlab version, no avi-movie will be saved')
    else
        saveavi=1;
        inmap = ls;
        xxx = strmatch(Wavename,inmap);
        if ~isempty(xxx)
            chosefrom='a':'z' ;
            Wavename = [chosefrom(randi(26,8,1)) '.avi'];
            disp(['Video exists, new name ' Wavename])
        end
    end
else
    saveavi=0;
end
saveavi=logical(saveavi);

figNo = gcf;
figure(figNo)  
clf
set(gca,'FontSize',14)
if nargin<2 || isempty(s)
  s=1;
end
Mv=[];
%disp('  Plotting frame by frame to record the movie...')
if ~ismatrix(Y.Z),
  [~,~,Nt]=size(Y.Z);
  if s==1
      axis equal
      colormap('winter');
      ax = [Y.y(1) Y.y(end) Y.x(1) Y.x(end) 3*min(Y.Z(:)) 3*max(Y.Z(:))];  
    for j=1:Nt
      set(gca,'nextplot','replacechildren');
      axis(ax)
      %axis('equal')
      surfl(Y.y,Y.x,2*Y.Z(:,:,j),[-30, 45]);
      shading interp
      xlabel('[m]')
      ylabel('[m]')
      view(-37.5,20)
      if isempty(Mv)
        Mv = getframe(figNo);
      else
        Mv(j)=getframe(figNo);
      end
    end  
    if saveavi,
        disp(['Saving ' Wavename]);
        savemovie(Mv,Wavename); %,'compression','none','fps',4*round(1/(Y.t(2)-Y.t(1))));
    end
        
  elseif s==2
    set(gca,'nextplot','replacechildren');
    box on
    for j=1:Nt
      contour(Y.y,Y.x,Y.Z(:,:,j),[0 0],'b')
      axis equal 
      xlabel('[m]')
      ylabel('[m]')
      if isempty(Mv)
        Mv = getframe(figNo);
      else
        Mv(j)=getframe(figNo);
      end  
    end
    if saveavi,
        disp(['Saving ' Wavename]);
        H = savemovie(Mv,Wavename); %,'compression','none','fps',4*round(1/(Y.t(2)-Y.t(1))));
    end
    
  else % s=3
    set(gca,'nextplot','replacechildren');
    box on
    colormap('gray')
    miz=min(Y.Z(:));
    maz=max(Y.Z(:));
    for j=1:Nt
      pcolor(Y.y,Y.x,Y.Z(:,:,j))
      caxis([miz maz])
      xlabel('[m]')
      ylabel('[m]')
      shading interp
      axis equal 
      axis([Y.y(1) Y.y(end) Y.x(1) Y.x(end)])% miz maz])
      if isempty(Mv)
        Mv = getframe(figNo);
      else
        Mv(j)=getframe(figNo);
      end  
    end
    hold on
    contour(Y.y,Y.x,Y.Z(:,:,end)+100,[100 100])
    hold off
    Mv(Nt+1:Nt+10)=getframe(figNo);
    if saveavi,
        H = savemovie(Mv,Wavename); %,'compression','none')
    end
  end    
elseif ndims(Y.Z)>1 && isfield(Y,'t')
  [~,Nt]=size(Y.Z);
  for j=1:Nt
    plot(Y.x,Y.Z(:,j))
    hold on
    plot([Y.x(1) Y.x(end)],[0 0],':')
    hold off
    xlabel('[m]')
    ylabel('[m]')
    axis([Y.x(1) Y.x(end),min(Y.Z(:))*2,max(Y.Z(:))*2])
    if isempty(Mv)
        Mv = getframe(figNo);
      else
        Mv(j)=getframe(figNo);
    end
  end
    if saveavi,
        H = savemovie(Mv,Wavename);% ,'compression','none')
    end
else
  if ~isfield(Y,'t')
    disp(...
'Can not make a movie without time variable, field .t must exist')
    return
  else
    disp('Wrong dimension of input. Can not make a movie')  
    return
  end
end
end

function HH = savemovie(ThisMovie,MovieName)
v = VideoWriter(MovieName);
open(v);
writeVideo(v,ThisMovie);
close(v)
HH = ['Movie' MovieName 'saved'];
end

