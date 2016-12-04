function varargout = kde1dgui(varargin)
%KDE1DGUI GUI to Kernel Density Estimator.
%
%  CALL:    kdegui;
%
%  Example
%   data = rndray(1,100,1);
%  if ismatlab,
%   kde1dgui
%  end
% 
% See also  kde2dgui, kde


% Tested on matlab7
% By pab Jan2005
% revised pab sept 2005
% - fixed some bugs + updated help header

% Edit the above text to modify the response to help kde1dgui

% Last Modified by GUIDE v2.5 19-Apr-2005 00:01:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @kde1dgui_OpeningFcn, ...
                   'gui_OutputFcn',  @kde1dgui_OutputFcn, ...
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


% --- Executes just before kde1dgui is made visible.
function kde1dgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for kde1dgui
handles.output = hObject;


handles.histPanel.bins = 10;

%handles.dataSelected = [];


kdeopts           = kdeoptset('kde');
kdeopts.fixh      = 1;
kdeopts.adaptive  = 0;
kdeopts.transform = 0;
kdeopts.addbumps  = 0;
kdeopts.bumpLinestyle = 'r--';
kdeopts.linestyle = 'k-';

%handles.kdePanel.kdeopts = kdeopts;

kdePanel.kdeopts = kdeopts;

setappdata(handles.figure1,'kdePanel',kdePanel)
% Populate kdePanel
updateKdePanel(handles);

setappdata(handles.figure1,'dataSelected',[])

% Populate the data menu
updateDataMenu(handles)
set(handles.dataMenu,'Value',1)
%Resample panel
set(handles.rbNo,'Value',1)


% Update handles structure
guidata(hObject, handles);


% UIWAIT makes kde1dgui wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function updateKdePanel(handles)
%UPDATEKDEPANEL Updates KDEPANEL of gui with current settings in kdePanel struct  
  

kdePanel = getappdata(handles.figure1,'kdePanel');
opts     = kdePanel.kdeopts;

%opts = handles.kdePanel.kdeopts;


if opts.adaptive==1 || opts.transform==1
   opts.addbumps=0;
   kdePanel.kdeopts.addbumps=0;
   set(handles.cbAddBumps,'Enable','off')
   set(handles.eBumpLinestyle,'Enable','off')
   % Update kdePanel structure
   setappdata(handles.figure1,'kdePanel',kdePanel);
   % Update handles structure
   %guidata(hObject, handles);
else
   set(handles.eBumpLinestyle,'Enable','on')
   set(handles.cbAddBumps,'Enable','on')
end
if kdePanel.kdeopts.addbumps ==1
  set(handles.eBumpLinestyle,'Visible','on')
else
  set(handles.eBumpLinestyle,'Visible','off')
end
if opts.fixh==1
  set(handles.eHs,'Enable','off')
  set(handles.hruleMenu,'Enable','on')
else
  set(handles.hruleMenu,'Enable','off')
  set(handles.eHs,'Enable','on')
end

set(handles.eHs,'String',num2str(opts.hs))
set(handles.eL2,'String',num2str(opts.L2))
set(handles.eAlpha,'String',num2str(opts.alpha))
set(handles.cbHRule, 'Value',opts.fixh);
set(handles.cbAdaptive, 'Value',opts.adaptive);
set(handles.cbTransform, 'Value',opts.transform)
set(handles.cbAddBumps, 'Value',opts.addbumps)
set(handles.eBumpLinestyle,'String',opts.bumpLinestyle)
set(handles.eLinestyle,'String',opts.linestyle)

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
function varargout = kde1dgui_OutputFcn(hObject, eventdata, handles)
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
   kdePanel = getappdata(handles.figure1,'kdePanel');
   kdePanel.kdeopts.hsMethod = strList{val};
   setappdata(handles.figure1,'kdePanel',kdePanel);
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
kdePanel = getappdata(handles.figure1,'kdePanel');
kdePanel.kdeopts.hs = hs;
setappdata(handles.figure1,'kdePanel',kdePanel);
   



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
kdePanel = getappdata(handles.figure1,'kdePanel');

kdePanel.kdeopts.transform=nonLinear;
kdePanel.kdeopts.L2 = L2;

setappdata(handles.figure1,'kdePanel',kdePanel)


updateKdePanel(handles)



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
kdePanel = getappdata(handles.figure1,'kdePanel');
kdePanel.kdeopts.adaptive=adaptive;
kdePanel.kdeopts.alpha = alpha1;

setappdata(handles.figure1,'kdePanel',kdePanel)
updateKdePanel(handles)



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
kdePanel = getappdata(handles.figure1,'kdePanel');
kdePanel.kdeopts.adaptive=adaptive;
kdePanel.kdeopts.alpha = alpha1;

setappdata(handles.figure1,'kdePanel',kdePanel);
updateKdePanel(handles)


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
  data = getappdata(handles.figure1,'dataSelected');
  if any(data<0)
    set(hObject,'Value',0);
    errordlg('Make sure data>0 for this option!','Error transform','modal')
    return
  end
   L2 = 0.5;
else
   L2=1;
end
nonLinear = (L2~=1);

kdePanel = getappdata(handles.figure1,'kdePanel');
kdePanel.kdeopts.transform=nonLinear;
kdePanel.kdeopts.L2 = L2;
setappdata(handles.figure1,'kdePanel',kdePanel)

updateKdePanel(handles)



% --- Executes on button press in cbHRule.
function cbHRule_Callback(hObject, eventdata, handles)
% hObject    handle to cbHRule (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbHRule

fixh = get(hObject,'Value');

if fixh==1
    set(handles.hruleMenu,'Enable','on')
   set(handles.eHs,'Enable','off')
else
   set(handles.hruleMenu,'Enable','off')
   set(handles.eHs,'Enable','on')
end
kdePanel = getappdata(handles.figure1,'kdePanel');

kdePanel.kdeopts.fixh=fixh;
setappdata(handles.figure1,'kdePanel',kdePanel);


% --- Executes on button press in cbAddBumps.
function cbAddBumps_Callback(hObject, eventdata, handles)
% hObject    handle to cbAddBumps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbAddBumps

kdePanel = getappdata(handles.figure1,'kdePanel');
kdePanel.kdeopts.addbumps=get(hObject,'Value');
if kdePanel.kdeopts.addbumps ==1
  set(handles.eBumpLinestyle,'Visible','on')
else
  set(handles.eBumpLinestyle,'Visible','off')
end
setappdata(handles.figure1,'kdePanel',kdePanel);
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

kdePanel = getappdata(handles.figure1,'kdePanel');
kdePanel.kdeopts.kernel = strList{val};
setappdata(handles.figure1,'kdePanel',kdePanel);

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
%UPDATEDATAMENU

% Updates the listbox to match the current workspace
vars = evalin('base','who');

for ix =length(vars):-1:1
  NsizeVar1 = sprintf('size(%s)',vars{ix});
  Nsize = evalin('base',NsizeVar1,'errordlg(lasterr,''Error generating plots'',''modal'')');
  isnumericVar1 = sprintf('isnumeric(%s)',vars{ix});
  isNumber = evalin('base',isnumericVar1,'errordlg(lasterr,''Error generating plots'',''modal'')');
  if prod(Nsize)~= max(Nsize) || ~isNumber 
    % remove data if dimension>1 or if not numeric data
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

kdePanel = getappdata(handles.figure1,'kdePanel');
buttonState = get(hObject,'Value') ;
if buttonState == get(hObject,'Max')
   
  ls = findNextLinestyle(kdePanel.kdeopts.linestyle);
  set(handles.tbHold,'String', 'Hold On');
  hold(handles.kdeAxes,'on')
   
elseif buttonState== get(hObject,'Min')
   
   ls = 'k-';
   set(handles.tbHold,'String', 'Hold Off');
   hold(handles.kdeAxes,'off')
end
set(handles.eLinestyle,'string',ls)
kdePanel.kdeopts.linestyle = ls;
setappdata(handles.figure1,'kdePanel',kdePanel);



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

data = getappdata(handles.figure1,'dataSelected');

if isempty(data)
  uiwait(msgbox('Load data before plotting!','No data loaded','modal'));
  return
end

kdePanel = getappdata(handles.figure1,'kdePanel');
if kdePanel.kdeopts.fixh==1
   kdePanel.kdeopts.hs = 0;
end

figure(handles.figure1)



pdf = kde(data,kdePanel.kdeopts);

pdf.title = sprintf('h = %s',num2str(pdf.options.hs));
%axis(handles.kdeAxes);
ls = kdePanel.kdeopts.linestyle;
if isempty(ls)
   pdfplot(pdf)
   %H=plot(pdf.x{1},pdf.f);
else
   pdfplot(pdf,ls)
   %H=plot(pdf.x{1},pdf.f,ls);
end

if ((kdePanel.kdeopts.addbumps==1) && ...
   (kdePanel.kdeopts.adaptive==0) && ...
   (kdePanel.kdeopts.transform==0))
   hs = pdf.options.hs;
   n  = length(data);
   x0 = linspace(-2,2);
   y  = mkernel(x0,pdf.options.kernel)/(n*hs); 
   hold_state = ishold;
   hold on
   linestyle = kdePanel.kdeopts.bumpLinestyle;
   for i=1:n
     plot(data(i)+x0*hs,y,linestyle);
     plot(repmat(data(i),1,2), [0 max(y)],linestyle)
   end
   if ~hold_state, 
      hold off,
   end
end

%pdfplot(pdf,ls)
hs = kdePanel.kdeopts.hs;
if ( isempty(hs) || hs<=0)
   kdePanel.kdeopts.hs=pdf.options.hs;
   setappdata(handles.figure1,'kdePanel',kdePanel)
   updateKdePanel(handles)
end

if ishold
  ls =   kdePanel.kdeopts.linestyle;
  nextLs = findNextLinestyle(ls);
  kdePanel.kdeopts.linestyle = nextLs;
  set(handles.eLinestyle,'String',nextLs)
  setappdata(handles.figure1,'kdePanel',kdePanel);
end

function nextLs = findNextLinestyle(ls)
  defaultColorOrder = 'kbrgmcy';
  defaultLinestyleOrder = {'-','--','-.',':',':.',':o',':x',':+',':*',':s'};
  nextColorIdx = 1;
  nextLsIdx    = 1;
  if ~isempty(ls)
    Nc             = length(defaultColorOrder); % Number of valid colors
    Nls            = length(defaultLinestyleOrder);
    startIndex     = double('a');
    endIndex       = double('z');
    shiftIndex     = startIndex-1;
    lettersHandled = endIndex-shiftIndex;

    ind1       = double(defaultColorOrder)- shiftIndex;
    ind2       = repmat(Nc+1,lettersHandled,1);
    ind2(ind1) = 1:Nc;

    ind3  = double(ls)-shiftIndex;

    ind = find((ind3<0|ind3>lettersHandled));
    if any(ind)
      ind3(ind) = lettersHandled;
    end
    ind = ind2(ind3);
    ix = find(ind==Nc+1); % indices to linestyle characters
    if any(ix)
      iz = strmatch(ls(ix),defaultLinestyleOrder,'exact');
      if any(iz)
        nextLsIdx =iz(1);
      end
    end
    iy = find(ind<=Nc); % index to color character
    if any(iy)
      nextColorIdx = ind(iy)+1;
      if nextColorIdx>Nc
        nextColorIdx = 1;
        nextLsIdx = mod(nextLsIdx,Nls)+1;
      end
    end
  end
  nextLs = [defaultColorOrder(nextColorIdx) defaultLinestyleOrder{nextLsIdx}];
    
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
      str = {' KDE1DGUI gives a GUI to compute the kernel density estimate',...
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

      

% --- Executes on key press over dataMenu with no controls selected.
function dataMenu_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to dataMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateDataMenu(handles)

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over dataMenu.
function dataMenu_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to dataMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateDataMenu(handles)

% --- Executes on selection change in dataMenu.
function dataMenu_Callback(hObject, eventdata, handles)
% hObject    handle to dataMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns dataMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dataMenu

updateDataMenu(handles)
loadData(handles)

% --- Executes on button press in pbLoadData.
function pbLoadData_Callback(hObject, eventdata, handles)
% hObject    handle to pbLoadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateDataMenu(handles)
loadData(handles)


function loadData(handles)
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
   dataSelected = sample(evalin('base',var1,'errordlg(lasterr,''Error generating plots'',''modal'')'),...
      Nsiz,doReplace);
else
   dataSelected = evalin('base',var1,'errordlg(lasterr,''Error generating plots'',''modal'')');
end
Nsize = size(dataSelected);
if prod(Nsize)~= max(Nsize)
   errordlg('Data selected is not a vector','Incorrect Selection','modal')
else
   set(handles.eSampleSize,'string',int2str(length(dataSelected)));
   setappdata(handles.figure1,'dataSelected',dataSelected);
   %guidata(hObject,handles)
end



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
function eBumpLinestyle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eBumpLinestyle (see GCBO)
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
   set(handles.rbReplace,'Value',0)
   dataSelected = getappdata(handles.figure1,'dataSelected');
   N = length(dataSelected);
   if N>0
     set(handles.eSampleSize,'String',int2str(N))
   end
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
   set(handles.rbReplace,'Value',0)
end

% --- Executes on button press in rbReplace.
function rbReplace_Callback(hObject, eventdata, handles)
% hObject    handle to rbReplace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbReplace
if get(hObject,'Value')==0
  dataSelected = getappdata(handles.figure1,'dataSelected');
  N =  length(dataSelected);
   if N>0,
      set(handles.eSampleSize,'String',int2str(N))
   end
else
  set(handles.rbNo,'Value',0)
  set(handles.rbYes,'Value',1)
end


function eLinestyle_Callback(hObject, eventdata, handles)
% hObject    handle to eLinestyle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eLinestyle as text
%        str2double(get(hObject,'String')) returns contents of eLinestyle as a double

kdePanel = getappdata(handles.figure1,'kdePanel');
kdePanel.kdeopts.linestyle = get(hObject,'String');

setappdata(handles.figure1,'kdePanel',kdePanel);


function eBumpLinestyle_Callback(hObject, eventdata, handles)
% hObject    handle to eBumpLinestyle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eBumpLinestyle as text
%        str2double(get(hObject,'String')) returns contents of eBumpLinestyle as a double

kdePanel = getappdata(handles.figure1,'kdePanel');
kdePanel.kdeopts.bumpLinestyle = get(hObject,'String');
setappdata(handles.figure1,'kdePanel',kdePanel);


% --- Executes on button press in plotHist.
function pbPlotHist_Callback(hObject, eventdata, handles)
% hObject    handle to plotHist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%val = get(hObject,'Value')

dataSelected = getappdata(handles.figure1,'dataSelected');

if ~isempty(dataSelected)
   figure(handles.figure1)
   N = handles.histPanel.bins;
   % newplot(handles.kdeAxes)
   histgrm(dataSelected,N,0,1)
end



