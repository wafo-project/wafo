function Mv=seamovie(Y,s)
% SEAMOVIE Makes a movie of a 2D (x,t) or 3D (x,y,t) simulated sea 
%  
% CALL:   Mv = seamovie(Y,s)
%  
%         Mv = movie
%         Y = struct with 2d or 3d simulation (from seasim)
%         s = type of plot if 3d: if s=1 then surf-plot, if s=2 contour,
%             else gray-scale overview with troughs dark and crests light
%             (default 1)
% 
% The recording is not very fast, each frame must be plotted and then saved
% using getframe. This may take a couple of seconds per frame.
% After the recording the movie is run, the number of frames per second is
% intended to be the true number as given by the time scale in input, but
% the resulting rate may vary depending on computer and network.
%  
% NB! Users with older Matlab than 5.3 may need to do some changes to the
% routine, see >> help getframe.  
%  
% Example: 
%   Y=seasim(demospec('dir'),2^8,1,20,10,[],.5,2);
%   Mv=seamovie(Y);
%  
% See also  seasim, movie, getframe
  
% Tested on Matlab 5.3 
% revised pab June 2005
% -fixed a bug: in matlab7: "Mv =[]; Mv(j) = getframe;" does not work, now
% fixed.
% revised es 20.06.00 if wrong dimension then message and return, not error
% Revised by es 13.06.00 more dimension checks  
% By es 23.05.00

figNo = gcf;
figure(figNo)  
if nargin<2||isempty(s)
  s=1;
end
Mv=[];
disp('  Plotting frame by frame to record the movie...')
if ndims(Y.Z)>2
  [Ny,Nx,Nt]=size(Y.Z);
  if s==1
    for j=1:Nt
      colormap('winter')
      surfl(Y.x,Y.y,Y.Z(:,:,j),[-30, 45]);
      shading interp
      view(-37.5,20)
      axis([Y.x(1) Y.x(end) Y.y(1) Y.y(end) 7*min(Y.Z(:)) 7*max(Y.Z(:))])
      set(gca,'xtick',[])   
      set(gca,'ytick',[])
      axis('square')
      axis('off')
      if isempty(Mv)
        Mv = getframe;
      else
        Mv(j)=getframe;
      end
    end
  elseif s==2
    for j=1:Nt
      contour(Y.x,Y.y,Y.Z(:,:,j),[0 0],'b')
      axis square
      xlabel('[m]')
      ylabel('[m]')
      if isempty(Mv)
        Mv = getframe(figNo);
      else
        Mv(j)=getframe(figNo);
      end
      
    end
  else
    colormap('gray')
    miz=min(Y.Z(:));
    maz=max(Y.Z(:));
    for j=1:Nt
      pcolor(Y.Z(:,:,j))
      caxis([miz maz])
      xlabel('[m]')
      ylabel('[m]')
      shading interp
      axis square %equal
      if isempty(Mv)
        Mv = getframe(figNo);
      else
        Mv(j)=getframe(figNo);
      end
      
    end
  end    
elseif ndims(Y.Z)>1 && isfield(Y,'t')
  [Nx,Nt]=size(Y.Z);
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
else
  if ~isfield(Y,'t')
    disp(...
'Can not make a movie without time variable, field .t must exist in input')
    return
  else
    disp('Wrong dimension of input. Can not make a movie')  
    return
  end
end
try
  disp('  Running the movie')
  clf
   movie(Mv,1,Y.t(2)-Y.t(1))
 catch
  for ix=1:length(Mv),
    t0 = clock;
    colormap(Mv(ix).colormap);
    image(Mv(ix).cdata),
    drawnow,shg,
    t1 = etime(clock,t0);
    ptime = Y.t(2)-Y.t(1) - t1;
    if ptime>0
      pause(ptime)
    end
  end
end
  
