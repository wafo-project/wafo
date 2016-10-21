function mkdeleteproxy(cax,ha,tag)
%MKDELETEPROXY  Make DeleteProxy object
%
%  CALL mkdeleteproxy(cax,ha,tag)
%
%   cax = handle to the axes the annotation objects should be bound to.
%   ha  = handle to annotation objects
%
% MKDELETEPROXY create an DeleteProxy object (an invisible text object in
% the axes, CAX) so that the other annotations objects will be deleted
% properly.

%History
% By pab March 2007

if nargin<3
  tag = 'DeleteProxy';
end

Nh = length(ha);
DeleteProxy = zeros(1,1+Nh);

DeleteProxy(1) = text('parent',cax,'visible','off',...
    'tag',tag,...
    'handlevisibility','off',...
    'deletefcn',@deletefun);
  
for ix = 1:Nh
  if strcmpi(get(ha,'type'),'text')
    DeleteProxy(ix+1) = ha(ix);
    set(ha(ix),'deletefcn',@deletefun)
  else
  DeleteProxy(ix+1) = text('parent',ha(ix),'visible','off',...
    'tag',tag,...
    'handlevisibility','off',...
    'deletefcn',@deletefun);
  end
end
  
set(DeleteProxy(1),'userdata',ha);
for ix = 2:Nh+1
  set(DeleteProxy(ix),'userdata',DeleteProxy([1:ix-1 ix+1:end]));
end

function deletefun(varargin)
% h = gcbo
for ix = 1:nargin
  delete(get(varargin{ix},'userdata'))
end
