function varargout = kde2dgui(varargin)
%KDE2DGUI GUI to Kernel Density Estimator in two dimensions.
%
%  CALL:  kde2dgui
%
% Example:
%
%  data = rndray(1,1000,2);
%  if ismatlab,
%    kde2dgui
%  end
%
% See also  kde1dgui, kde

% Tested on matlab 7.0
% By pab June2005
%  Revised pab sept 2005
% fixed some bugs + updated help header

% Edit the above text to modify the response to help kde2dgui

% Last Modified by GUIDE v2.5 05-Mar-2005 22:08:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @kde2dgui_OpeningFcn, ...
                   'gui_OutputFcn',  @kde2dgui_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before kde2dgui is made visible.
function kde2dgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for kde2dgui
handles.output = hObject;
kdeopts           = kdeoptset('kde');
kdeopts.fixh      = 1;
kdeopts.adaptive  = 0;
kdeopts.transform = 0;
%kdeopts.addbumps  = 0;
%kdeopts.bumpLinestyle = 'r--';
%kdeopts.linestyle = 'k-';

handles.kdePanel.kdeopts = kdeopts;


handles.kdefigure = figure;
set(handles.kdefigure,'Tag','kde2dgui')
handles.kdeAxes = axes;


handles.dataSelected = [];
% Update handles structure
guidata(hObject, handles);

%Resample panel
set(handles.rbNo,'Value',1)

% Populate the data menu
updateDataMenu(handles)
set(handles.dataMenu,'Value',1)
% Populate kdePanel
updateKdePanel(handles,hObject);


% UIWAIT makes kde2dgui wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function updateKdePanel(handles,hObject)
opts = handles.kdePanel.kdeopts;

if opts.fixh==1
   set(handles.eHs,'Enable','off')
   set(handles.hruleMenu,'Enable','on')
else
    set(handles.hruleMenu,'Enable','off')
   set(handles.eHs,'Enable','on')
end

set(handles.eHs,'String',sprintf('[%s]',num2str(opts.hs)))
set(handles.eL2,'String',sprintf('[%s]',num2str(opts.L2)))
set(handles.eAlpha,'String',sprintf('[%s]',num2str(opts.alpha)))
set(handles.cbHRule, 'Value',opts.fixh);
set(handles.cbAdaptive, 'Value',opts.adaptive);
set(handles.cbTransform, 'Value',opts.transform)
%set(handles.cbAddScatter, 'Value',opts.addbumps)
%set(handles.eScatterLinestyle,'String',opts.bumpLinestyle)
%set(handles.eLinestyle,'String',opts.linestyle)

kernelEntries = lower(get(handles.kernelMenu,'String'));
ix = strmatch(lower(opts.kernel),kernelEntries);
if length(ix) == 1
    set(handles.kernelMenu,'Value',ix);
end 
hruleEntries = get(handles.hruleMenu,'String');
ix = strmatch(opts.hsMethod,hruleEntries);
if length(ix) == 1
    set(handles.hruleMenu,'Value',ix);
end 


% --- Outputs from this function are returned to the command line.
function varargout = kde2dgui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes on selection change in hruleMenu.
function hruleMenu_Callback(hObject, eventdata, handles)
% hObject    handle to hruleMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns hruleMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from hruleMenu

%if get(handles.cbHRule,'Value')==0
   val     = get(hObject,'Value');
   strList = get(hObject,'String');
   handles.kdePanel.kdeopts.hsMethod = strList{val};
   guidata(hObject, handles);
%end


% --- Executes during object creation, after setting all properties.
function hruleMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hruleMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function eHs_Callback(hObject, eventdata, handles)
% hObject    handle to eHs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eHs as text
%        str2double(get(hObject,'String')) returns contents of eHs as a double

hs = str2double(get(hObject,'String'));

if (isnan(hs))
   hs = [];
end
handles.kdePanel.kdeopts.hs = hs;
guidata(hObject,handles)
   



% --- Executes during object creation, after setting all properties.
function eHs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eHs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function eL2_Callback(hObject, eventdata, handles)
% hObject    handle to eL2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eL2 as text
%        str2double(get(hObject,'String')) returns contents of eL2 as a double
L2 = str2double(get(hObject,'String'));
invalidInput = isnan(L2)| (L2 <0);
if (invalidInput) 
   L2 = 1;
end
nonLinear = (L2~=1);
handles.kdePanel.kdeopts.transform=nonLinear;

handles.kdePanel.kdeopts.L2 = L2;
guidata(hObject,handles)

updateKdePanel(handles,hObject)



% --- Executes during object creation, after setting all properties.
function eL2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eL2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cbAdaptive.
function cbAdaptive_Callback(hObject, eventdata, handles)
% hObject    handle to cbAdaptive (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbAdaptive

val = get(hObject,'Value');
if val==1
   alpha1 = 0.5;
else
   alpha1=0;
end
adaptive = (alpha1~=0);
handles.kdePanel.kdeopts.adaptive=adaptive;

handles.kdePanel.kdeopts.alpha = alpha1;
guidata(hObject,handles)

updateKdePanel(handles,hObject)



function eAlpha_Callback(hObject, eventdata, handles)
% hObject    handle to eAlpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eAlpha as text
%        str2double(get(hObject,'String')) returns contents of eAlpha as a double

alpha1 = str2double(get(hObject,'String'));
invalidInput = isnan(alpha1) | (alpha1<0) | (1<alpha1);
if (invalidInput) 
   alpha1 = 0;
end
adaptive = (alpha1~=0);
handles.kdePanel.kdeopts.adaptive=adaptive;

handles.kdePanel.kdeopts.alpha = alpha1;
guidata(hObject,handles)
updateKdePanel(handles,hObject)


% --- Executes during object creation, after setting all properties.
function eAlpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eAlpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cbTransform.
function cbTransform_Callback(hObject, eventdata, handles)
% hObject    handle to cbTransform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbTransform


val = get(hObject,'Value');
if val==1
   L2 = 0.5;
else
   L2=1;
end
nonLinear = (L2~=1);
handles.kdePanel.kdeopts.transform=nonLinear;

handles.kdePanel.kdeopts.L2 = L2;
guidata(hObject,handles)

updateKdePanel(handles,hObject)



% --- Executes on button press in cbHRule.
function cbHRule_Callback(hObject, eventdata, handles)
% hObject    handle to cbHRule (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbHRule

fixh = get(hObject,'Value');
handles.kdePanel.kdeopts.fixh=fixh;
if fixh==1
    set(handles.hruleMenu,'Enable','on')
   set(handles.eHs,'Enable','off')
else
   set(handles.hruleMenu,'Enable','off')
   set(handles.eHs,'Enable','on')
end
guidata(hObject,handles)

%updateKdePanel(handles)


% --- Executes on button press in cbAddScatter.
function cbAddScatter_Callback(hObject, eventdata, handles)
% hObject    handle to cbAddScatter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbAddScatter

handles.kdePanel.kdeopts.addbumps=get(hObject,'Value');
guidata(hObject,handles)
%updateKdePanel(handles)

% --- Executes on selection change in kernelMenu.
function kernelMenu_Callback(hObject, eventdata, handles)
% hObject    handle to kernelMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns kernelMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from kernelMenu

val     = get(hObject,'Value');
strList = get(hObject,'String');
handles.kdePanel.kdeopts.kernel = strList{val};
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function kernelMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kernelMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function eNumBins_Callback(hObject, eventdata, handles)
% hObject    handle to eNumBins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eNumBins as text
%        str2double(get(hObject,'String')) returns contents of eNumBins as a double

bins = str2double(get(hObject,'String'));
invalidInput = isnan(bins) | (bins<1);
if (invalidInput) 
   bins = 10;
   set(hObject,'String',int2str(bins))
end

handles.histPanel.bins=bins;

guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function eNumBins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eNumBins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function eSampleSize_Callback(hObject, eventdata, handles)
% hObject    handle to eSampleSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eSampleSize as text
%        str2double(get(hObject,'String')) returns contents of eSampleSize as a double


% --- Executes during object creation, after setting all properties.
function eSampleSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eSampleSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function updateDataMenu(handles)
% hObject    handle to update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1

% Updates the listbox to match the current workspace
vars = evalin('base','who');
for ix =length(vars):-1:1
  NsizeVar1 = sprintf('size(%s)',vars{ix});
  Nsize = evalin('base',NsizeVar1,'errordlg(lasterr,''Error generating plots'',''modal'')');
  isnumericVar1 = sprintf('isnumeric(%s)',vars{ix});
  isNumber = evalin('base',isnumericVar1,'errordlg(lasterr,''Error generating plots'',''modal'')');
  if length(Nsize)>2 || Nsize(2)~= 2 ||  ~isNumber 
    % remove data if dimension~=2 or if not numeric data
    vars(ix) = [];
  end
end
if ~isempty(vars)
  set(handles.dataMenu,'String',vars)
end


% --- Executes on button press in tbHold.
function tbHold_Callback(hObject, eventdata, handles)
% hObject    handle to tbHold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tbHold

buttonState = get(hObject,'Value') ;
if buttonState == get(hObject,'Max')
   set(handles.tbHold,'String', 'Hold On');
   hold(handles.kdeAxes,'on')
elseif buttonState== get(hObject,'Min')
   set(handles.tbHold,'String', 'Hold Off');
   hold(handles.kdeAxes,'off')
end




% --- Executes on button press in pbClose.
function pbClose_Callback(hObject, eventdata, handles)
% hObject    handle to pbClose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(gcbf)


% --- Executes on button press in pbPlotKde.
function pbPlotKde_Callback(hObject, eventdata, handles)
% hObject    handle to pbPlotKde (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%data = evalin('base',handles.lbDataSelected,'errordlg(lasterr,''Error generating plots'',''modal'')');

if isempty(handles.dataSelected)
  uiwait(msgbox('Load data before plotting!','No data loaded','modal'));
   return
end
if handles.kdePanel.kdeopts.fixh==1
   handles.kdePanel.kdeopts.hs = 0;
end

figs = findobj('Type', 'figure');
makeNewFigure = isempty(figs);
if ~makeNewFigure
   h = findobj(figs, 'Tag', 'kde2dgui');
   if isempty(h)
      h = figure;
      handles.kdeAxes = axes;
      guidata(hObject,handles)
   end
   if length(h)>1
     h = h(1); % chose the first figure
   end
else
   h = figure;
   handles.kdeAxes = axes;
   guidata(hObject,handles)
end
figure(h)
set(h,'Tag','kde2dgui')

%figure(handles.figure1)

pdf = kde(handles.dataSelected,handles.kdePanel.kdeopts);

pdf.title = sprintf('h = %s',num2str(pdf.options.hs));


 PL = str2num(get(handles.eContourLevels,'String'));
 if isempty(PL)
   if isfield(pdf,'pl') 
     set(handles.eContourLevels,'String',sprintf('%g ',pdf.pl))
   end
 else
   [ql, PL] = qlevels(pdf.f,PL);
   pdf.cl   = ql;
   pdf.pl   = PL;
 end


%axis(handles.kdeAxes);

listEntries = get(handles.plotMenu,'String');
plotflag = get(handles.plotMenu,'Value');
if length(plotflag) ~= 1
    errordlg('You must select one plot option','Incorrect Selection','modal')
end 
plotMetod = listEntries{plotflag};

pShading = get(handles.cbShading,'Value');
if pShading
  shadingtxt = 'interp'
else
  shadingtxt = 'faceted';
end

ls = get(handles.eLinestyle,'String');
if isempty(ls)
   pdfplot(pdf,plotflag,shadingtxt)
else
   pdfplot(pdf,ls,plotflag,shadingtxt)
end

plotScatter = get(handles.cbAddScatter,'Value');
if plotScatter
   hold_state = ishold;
   hold on
   linestyle = get(handles.eScatterLinestyle,'String');
   plot(handles.dataSelected(:,1),handles.dataSelected(:,2),linestyle)
   
   if ~hold_state, 
      hold off,
   end
end

%pdfplot(pdf,ls)
hs = handles.kdePanel.kdeopts.hs;
if ( isempty(hs) || any(hs<=0))
   handles.kdePanel.kdeopts.hs=pdf.options.hs;
   guidata(hObject,handles)
   updateKdePanel(handles,hObject)
end

% --- Executes on key press over pbClose with no controls selected.
function pbClose_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to pbClose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch(eventdata.Key)
 case {'return'}
  if ~strcmp(get(obj,'UserData'),'Cancel')
      set(gcbf,'UserData','OK');
      uiresume(gcbf);
  else
      delete(gcbf)
  end
 case 'escape'
  delete(gcbf)
end



% --- Executes on button press in pbHelpButton.
function pbHelpButton_Callback(hObject, eventdata, handles)
% hObject    handle to pbHelpButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

helpdlg(kdehelpstr,'kde1Dgui help')

   function str = kdehelpstr
      str = {' KDE2DGUI gives a GUI to compute the kernel density estimate',...
'  Notice that densities close to normality appear to be the easiest for the kernel', ...
'  estimator to estimate and that the degree of estimation difficulty increases with ', ...
'   skewness, kurtosis and multimodality.',...
'   If L2~=1 KDE transforms the data before estimation. The final estimate',...
'   is obtained by transforming back by a simple change of variables.',...
'   Beaware of spurious spikes close to the edges when L2~=1.',...
'   These spikes are due to numerical problems close to the edges.',...
' ',...
'    kernel = String defining the kernel function.',....
'    hs     = smooting parameter vector/matrix.',...
'             (default compute from data using hsMethod)',...
'  hsMethod = string defining the method to compute the  smooting',...
'             parameter hs',....
'    alpha  = sensitivity parameter ',...
'             A good choice might be alpha = 0.5 ( or 1/D)',...
'             alpha = 0      Regular  KDE (hs is constant)',...
'             0 < alpha <= 1 Adaptive KDE (Make hs change adaptively)  ',...
'    L2     = transformation parameter (L2=1 means no transformation)',...
'             t(xi;L2) = xi^L2*sign(L2)   for L2(i) ~= 0',...
'             t(xi;L2) = log(xi)          for L2(i) == 0 '};

      

% --- Executes on selection change in dataMenu.
function dataMenu_Callback(hObject, eventdata, handles)
% hObject    handle to dataMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns dataMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dataMenu

updateDataMenu(handles)
loadData(hObject,handles)

% --- Executes on button press in pbLoadData.
function pbLoadData_Callback(hObject, eventdata, handles)
% hObject    handle to pbLoadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateDataMenu(handles)
loadData(hObject,handles)


function loadData(hObject,handles)
%LOADDATA from workspace, possibly resampled
   listEntries = get(handles.dataMenu,'String');
index_selected = get(handles.dataMenu,'Value');
if length(index_selected) ~= 1
    errordlg('You must select one variable','Incorrect Selection','modal')
else
    var1 = listEntries{index_selected};
end 
doResample = get(handles.rbYes,'Value');

if (doResample==1) 
   doReplace = get(handles.rbReplace,'Value');
   Nsiz = str2double(get(handles.eSampleSize,'string'));
   handles.dataSelected = sample(evalin('base',var1,'errordlg(lasterr,''Error generating plots'',''modal'')'),...
      Nsiz,doReplace);
else
   handles.dataSelected = evalin('base',var1,'errordlg(lasterr,''Error generating plots'',''modal'')');
end
set(handles.eSampleSize,'string',int2str(length(handles.dataSelected)));
guidata(hObject,handles)




% --- Executes during object creation, after setting all properties.
function dataMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function eScatterLinestyle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eScatterLinestyle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function eLinestyle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eLinestyle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rbNo.
function rbNo_Callback(hObject, eventdata, handles)
% hObject    handle to rbNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbNo

val = get(hObject,'Value');
if val==1
   set(handles.rbYes,'Value',0)
else
   set(handles.rbYes,'Value',1)
end
%guidata(hObject,handles)

% --- Executes on button press in rbYes.
function rbYes_Callback(hObject, eventdata, handles)
% hObject    handle to rbYes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbYes

val = get(hObject,'Value');
if val==1
   set(handles.rbNo,'Value',0)
else
   set(handles.rbNo,'Value',1)
end

% --- Executes on button press in rbReplace.
function rbReplace_Callback(hObject, eventdata, handles)
% hObject    handle to rbReplace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbReplace
if get(hObject,'Value')==0
   if ~isempty(handles.dataSelected)
      N = length(handles.dataSelected)
      set(handles.eSampleSize,'String',int2str(N))
   end
end


function eLinestyle_Callback(hObject, eventdata, handles)
% hObject    handle to eLinestyle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eLinestyle as text
%        str2double(get(hObject,'String')) returns contents of eLinestyle as a double

handles.kdePanel.kdeopts.linestyle = get(hObject,'String');
%updateKdePanel(handles)
guidata(hObject,handles)



function eScatterLinestyle_Callback(hObject, eventdata, handles)
% hObject    handle to eScatterLinestyle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eScatterLinestyle as text
%        str2double(get(hObject,'String')) returns contents of eScatterLinestyle as a double

handles.kdePanel.kdeopts.bumpLinestyle = get(hObject,'String');
%updateKdePanel(handles)
guidata(hObject,handles)



% --- Executes on selection change in plotMenu.
function plotMenu_Callback(hObject, eventdata, handles)
% hObject    handle to plotMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns plotMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plotMenu


% --- Executes during object creation, after setting all properties.
function plotMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plotMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function eContourLevels_Callback(hObject, eventdata, handles)
% hObject    handle to eContourLevels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eContourLevels as text
%        str2double(get(hObject,'String')) returns contents of eContourLevels as a double


% --- Executes during object creation, after setting all properties.
function eContourLevels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eContourLevels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb2Dview.
function pb2Dview_Callback(hObject, eventdata, handles)
% hObject    handle to pb2Dview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figs = findobj('Type', 'figure');
if ~isempty(figs)
   h = findobj(figs, 'Tag', 'kde2dgui');
   if ~isempty(h)
     figure(h(1))
     view(2)
   end
end

% --- Executes on button press in pb3Dview.
function pb3Dview_Callback(hObject, eventdata, handles)
% hObject    handle to pb3Dview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figs = findobj('Type', 'figure');
if ~isempty(figs)
   h = findobj(figs, 'Tag', 'kde2dgui');
   if ~isempty(h)
     figure(h(1))
     view(3)
   end
end


% --- Executes on button press in cbShading.
function cbShading_Callback(hObject, eventdata, handles)
% hObject    handle to cbShading (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbShading


