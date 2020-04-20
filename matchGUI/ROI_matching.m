function varargout = ROI_matching(varargin)
%ROI_MATCHING MATLAB code file for ROI_matching.fig
%      ROI_MATCHING, by itself, creates a new ROI_MATCHING or raises the existing
%      singleton*.
%
%      H = ROI_MATCHING returns the handle to a new ROI_MATCHING or the handle to
%      the existing singleton*.
%
%      ROI_MATCHING('Property','Value',...) creates a new ROI_MATCHING using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to ROI_matching_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      ROI_MATCHING('CALLBACK') and ROI_MATCHING('CALLBACK',hObject,...) call the
%      local function named CALLBACK in ROI_MATCHING.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ROI_matching

% Last Modified by GUIDE v2.5 31-Aug-2018 16:47:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ROI_matching_OpeningFcn, ...
                   'gui_OutputFcn',  @ROI_matching_OutputFcn, ...
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


% --- Executes just before ROI_matching is made visible.
function ROI_matching_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for ROI_matching
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ROI_matching wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ROI_matching_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function entry_ROI_adjacency_Callback(hObject, eventdata, handles)
% hObject    handle to entry_ROI_adjacency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of entry_ROI_adjacency as text
%        str2double(get(hObject,'String')) returns contents of entry_ROI_adjacency as a double


% --- Executes on selection change in dropdown_filter_ll_gg.
function dropdown_filter_ll_gg_Callback(hObject, eventdata, handles)
% hObject    handle to dropdown_filter_ll_gg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns dropdown_filter_ll_gg contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dropdown_filter_ll_gg


% --- Executes during object creation, after setting all properties.
function dropdown_filter_ll_gg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dropdown_filter_ll_gg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function entry_filter_value_Callback(hObject, eventdata, handles)
% hObject    handle to entry_filter_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of entry_filter_value as text
%        str2double(get(hObject,'String')) returns contents of entry_filter_value as a double


% --- Executes on selection change in dropdown_filter_type.
function dropdown_filter_type_Callback(hObject, eventdata, handles)
% hObject    handle to dropdown_filter_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns dropdown_filter_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dropdown_filter_type


% --- Executes during object creation, after setting all properties.
function dropdown_filter_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dropdown_filter_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function entry_cluster_ID_Callback(hObject, eventdata, handles)
% hObject    handle to entry_cluster_ID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of entry_cluster_ID as text
%        str2double(get(hObject,'String')) returns contents of entry_cluster_ID as a double


% --- Executes on button press in button_prev_cluster_ID.
function button_prev_cluster_ID_Callback(hObject, eventdata, handles)
% hObject    handle to button_prev_cluster_ID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in button_next_cluster_ID.
function button_next_cluster_ID_Callback(hObject, eventdata, handles)
% hObject    handle to button_next_cluster_ID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in button_clustering.
function button_clustering_Callback(hObject, eventdata, handles)
% hObject    handle to button_clustering (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox_filter_active.
function checkbox_filter_active_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_filter_active (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_filter_active


% --- Executes on selection change in dropdown_filter2_ll_gg.
function dropdown_filter2_ll_gg_Callback(hObject, eventdata, handles)
% hObject    handle to dropdown_filter2_ll_gg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns dropdown_filter2_ll_gg contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dropdown_filter2_ll_gg


% --- Executes during object creation, after setting all properties.
function dropdown_filter2_ll_gg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dropdown_filter2_ll_gg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function entry_filter2_value_Callback(hObject, eventdata, handles)
% hObject    handle to entry_filter2_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of entry_filter2_value as text
%        str2double(get(hObject,'String')) returns contents of entry_filter2_value as a double


% --- Executes on selection change in dropdown_filter2_type.
function dropdown_filter2_type_Callback(hObject, eventdata, handles)
% hObject    handle to dropdown_filter2_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns dropdown_filter2_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dropdown_filter2_type


% --- Executes during object creation, after setting all properties.
function dropdown_filter2_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dropdown_filter2_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function entry_display_session_Callback(hObject, eventdata, handles)
% hObject    handle to entry_display_session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of entry_display_session as text
%        str2double(get(hObject,'String')) returns contents of entry_display_session as a double


% --- Executes on button press in button_prev_session.
function button_prev_session_Callback(hObject, eventdata, handles)
% hObject    handle to button_prev_session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in button_next_session.
function button_next_session_Callback(hObject, eventdata, handles)
% hObject    handle to button_next_session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox_show_all_sessions.
function checkbox_show_all_sessions_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_show_all_sessions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_show_all_sessions


% --- Executes on button press in checkbox_rotate3d.
function checkbox_rotate3d_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_rotate3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_rotate3d


% --- Executes on button press in checkbox_finished.
function checkbox_finished_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_finished (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_finished


% --- Executes on button press in button_multi_remove_IDs.
function button_multi_remove_IDs_Callback(hObject, eventdata, handles)
% hObject    handle to button_multi_remove_IDs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in button_discard_cluster.
function button_discard_cluster_Callback(hObject, eventdata, handles)
% hObject    handle to button_discard_cluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in button_single_remove_IDs.
function button_single_remove_IDs_Callback(hObject, eventdata, handles)
% hObject    handle to button_single_remove_IDs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in radio_split.
function radio_split_Callback(hObject, eventdata, handles)
% hObject    handle to radio_split (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_split


% --- Executes on button press in radio_merge.
function radio_merge_Callback(hObject, eventdata, handles)
% hObject    handle to radio_merge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_merge


% --- Executes on button press in radio_belong.
function radio_belong_Callback(hObject, eventdata, handles)
% hObject    handle to radio_belong (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_belong


% --- Executes on button press in radio_clusterdisplay_3D.
function radio_clusterdisplay_3D_Callback(hObject, eventdata, handles)
% hObject    handle to radio_clusterdisplay_3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_clusterdisplay_3D


% --- Executes on button press in radio_clusterdisplay_2D.
function radio_clusterdisplay_2D_Callback(hObject, eventdata, handles)
% hObject    handle to radio_clusterdisplay_2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_clusterdisplay_2D


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in button_get_single_ROI.
function button_get_single_ROI_Callback(hObject, eventdata, handles)
% hObject    handle to button_get_single_ROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
