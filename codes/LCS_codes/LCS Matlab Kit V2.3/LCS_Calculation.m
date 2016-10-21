function varargout = LCS_Calculation(varargin)
% LCS_Calculation M-file for LCS_Calculation.fig
%      LCS_Calculation, by itself, creates a new LCS_Calculation or raises the existing
%      singleton*.
%
%      H = LCS_Calculation returns the handle to a new LCS_Calculation or
%      the handle to
%      the existing singleton*.
%
%      LCS_Calculation('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LCS_Calculation.M with the given input arguments.
%
%      LCS_Calculation('Property','Value',...) creates a new LCS_Calculation or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LCS_Calculation_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LCS_Calculation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LCS_Calculation

% Last Modified by GUIDE v2.5 07-Mar-2009 21:59:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LCS_Calculation_OpeningFcn, ...
                   'gui_OutputFcn',  @LCS_Calculation_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before LCS_Calculation is made visible.
function LCS_Calculation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LCS_Calculation (see VARARGIN)

% Choose default command line output for LCS_Calculation
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes LCS_Calculation wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = LCS_Calculation_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function starting_frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to starting_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function starting_frame_Callback(hObject, eventdata, handles)
% hObject    handle to starting_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of starting_frame as text
%        str2double(get(hObject,'String')) returns contents of starting_frame as a double


% --- Executes during object creation, after setting all properties.
function end_frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to starting_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function end_frame_Callback(hObject, eventdata, handles)
% hObject    handle to end_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of end_frame as text
%        str2double(get(hObject,'String')) returns contents of end_frame as a
%        double


% --- Executes during object creation, after setting all properties.
function interval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function interval_Callback(hObject, eventdata, handles)
% hObject    handle to interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of interval as text
%        str2double(get(hObject,'String')) returns contents of interval as a double


% --- Executes during object creation, after setting all properties.
function integration_length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to integration_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function integration_length_Callback(hObject, eventdata, handles)
% hObject    handle to integration_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of integration_length as text
%        str2double(get(hObject,'String')) returns contents of integration_length as a double


% --- Executes during object creation, after setting all properties.
function time_interval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function time_interval_Callback(hObject, eventdata, handles)
% hObject    handle to time_interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time_interval as text
%        str2double(get(hObject,'String')) returns contents of time_interval as a double


% --- Executes during object creation, after setting all properties.
function meshsize_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to meshsize_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function meshsize_x_Callback(hObject, eventdata, handles)
% hObject    handle to meshsize_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of meshsize_x as text
%        str2double(get(hObject,'String')) returns contents of meshsize_x as a double


% --- Executes during object creation, after setting all properties.
function meshsize_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to meshsize_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function meshsize_y_Callback(hObject, eventdata, handles)
% hObject    handle to meshsize_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of meshsize_y as text
%        str2double(get(hObject,'String')) returns contents of meshsize_y as a double


% --- Executes during object creation, after setting all properties.
function region_x_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to region_x_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function region_x_min_Callback(hObject, eventdata, handles)
% hObject    handle to region_x_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of region_x_min as text
%        str2double(get(hObject,'String')) returns contents of region_x_min as a double


% --- Executes during object creation, after setting all properties.
function region_x_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to region_x_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function region_x_max_Callback(hObject, eventdata, handles)
% hObject    handle to region_x_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of region_x_max as text
%        str2double(get(hObject,'String')) returns contents of region_x_max as a double



% --- Executes during object creation, after setting all properties.
function region_y_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to region_y_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function region_y_max_Callback(hObject, eventdata, handles)
% hObject    handle to region_y_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of region_y_max as text
%        str2double(get(hObject,'String')) returns contents of region_y_max as a double


% --- Executes during object creation, after setting all properties.
function region_y_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to region_y_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function region_y_min_Callback(hObject, eventdata, handles)
% hObject    handle to region_y_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of region_y_min as text
%        str2double(get(hObject,'String')) returns contents of region_y_min as a double


% --- Executes during object creation, after setting all properties.
function step_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to step_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function step_size_Callback(hObject, eventdata, handles)
% hObject    handle to step_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of step_size as text
%        str2double(get(hObject,'String')) returns contents of step_size as a double


% --- Executes on button press in forward_FTLE.
function forward_FTLE_Callback(hObject, eventdata, handles)
% hObject    handle to forward_FTLE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of forward_FTLE


% --- Executes on button press in backward_FTLE.
function backward_FTLE_Callback(hObject, eventdata, handles)
% hObject    handle to backward_FTLE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of backward_FTLE

% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Add the current directory to the path, as the pwd might change thru' the
% gui. Remove the directory from the path when gui is closed 
% (See figure1_DeleteFcn)
setappdata(hObject, 'StartPath', pwd);
addpath(pwd);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Menu Control%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Menu Control.1.File%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Remove the directory added to the path in the figure1_CreateFcn.
if isappdata(hObject, 'StartPath')
    rmpath(getappdata(hObject, 'StartPath'));
end


% --- Executes during object creation, after setting all properties.
function FTLE_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to calculation_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object deletion, before destroying properties.
function FTLE_type_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to FTLE_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)

% --------------------------------------------------------------------
function NewMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to NewMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function SaveMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to SaveMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function SaveAsMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to SaveAsMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Menu Control.2.Input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --------------------------------------------------------------------
function InputMenu_Callback(hObject, eventdata, handles)
% hObject    handle to InputMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function LoadOriginalDataMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to LoadOriginalDataMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global originaldata_filename originaldata_pathname find_txt find_wk1
[originaldata_filename, originaldata_pathname, filterindex] = uigetfile({'*.txt';'*.wk1'}, 'MultiSelect', 'on');

if iscell(originaldata_filename)==1
    find_txt=max(size(strfind(originaldata_filename{1},'.txt')));
    find_wk1=max(size(strfind(originaldata_filename{1},'.wk1')));
else
    find_txt=max(size(strfind(originaldata_filename,'.txt')));
    find_wk1=max(size(strfind(originaldata_filename,'.wk1')));
end

if find_txt==0 & find_wk1==0
    warndlg('Selection failed. Unknown file format.','error');
end

if filterindex==0
    warndlg('Selection failed. Matlab sometimes fails to a open large number of files. Please select less files.','error');
end
originaldata_filename=sort_char_array(originaldata_filename);

% --------------------------------------------------------------------
function ConvertToMATFilesMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to ConvertToMATFilesMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global originaldata_filename originaldata_pathname directory_name find_txt find_wk1

directory_name = uigetdir('','Select directory of the folder to which MAT files are to be saved');

prompt = {'Enter the starting frame number:'};
dlg_title = 'Convert Original Data to Matlab Files';
num_lines = 1;
def = {'1'};

MatFileStartingFrameString = inputdlg(prompt,dlg_title,num_lines,def);
MatFileStartingFrame = str2num(MatFileStartingFrameString{1});

n_originaldata=length(originaldata_filename);
for i=1:n_originaldata
    i_frame=i+MatFileStartingFrame-1;
    filename=[originaldata_pathname originaldata_filename{i}];
    
    if find_txt==1
        A = dlmread(filename);
    end
    if find_wk1==1
        A = wk1read(filename);
    end
    A = sortrows(A,1);
    B=sortrows(A,2);
    ndata=max(size(A));

    if i==1
        ysize=max(find(A(:,1)==A(1,1)));
        xsize=max(find(B(:,2)==B(1,2)));
        xcoor=B(1:xsize,1);
        ycoor=A(1:ysize,2);
        [xloc yloc]=meshgrid(xcoor,ycoor);
        save([directory_name '/xgrid'],'xcoor');
        save([directory_name '/ygrid'],'ycoor');
    end

    for j=1:xsize
        u(1:ysize,j)=A(1+ysize*(j-1):ysize*j,3);
        v(1:ysize,j)=A(1+ysize*(j-1):ysize*j,4);
    end

    expression=['u' num2str(i_frame) '=u'';'];
    eval(expression);
    expression=['v' num2str(i_frame) '=v'';'];
    eval(expression);
    expression=['save(''' directory_name '/U_T' num2str(i_frame) ''',''u' num2str(i_frame) ''');'];
    eval(expression);
    expression=['save(''' directory_name '/V_T' num2str(i_frame) ''',''v' num2str(i_frame) ''');'];
    eval(expression);
    expression=['clear u' num2str(i_frame) ' v' num2str(i_frame) ';'];
    eval(expression);
end

% --------------------------------------------------------------------
function PlotVelocityFieldMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PlotVelocityFieldMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global directory_name

prompt = {'Enter starting frame:','Enter ending frame:'};
dlg_title = 'Enter frame number to plot velocity';
num_lines = 1;
def = {'40','60'};
answer = inputdlg(prompt,dlg_title,num_lines,def);

n_velocity_start=str2num(answer{1});
n_velocity_end=str2num(answer{2});

velocity_plot(directory_name,n_velocity_start,n_velocity_end,handles.axes1);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Menu Control.3.Calculation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------------------
function FTLECalculationMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FTLECalculationMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function LoadVelocityFilesMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to LoadVelocityFilesMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global directory_name
directory_name = uigetdir('','Select directory of the folder for velocity files');

% --------------------------------------------------------------------
function ForwardBackwardCalculationMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to ForwardBackwardCalculationMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
helpdlg('Please choose whether forward or backward calculation to be calculated.','Forward/Backward Calculation')


% --------------------------------------------------------------------
function FTLE_Frames_Callback(hObject, eventdata, handles)
% hObject    handle to Calculation_Frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
helpdlg('Please fill in the FTLE Frames edit boxes for values of calculation frames.',...
    'Calculation_Frames')

% --------------------------------------------------------------------
function IntegrationManagementMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to TimeManagementMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
helpdlg('Please fill in the Integration edit boxes for values of integration parameters.',...
    'Integration Management')


% --------------------------------------------------------------------
function CalculationRegionMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CalculationRegionMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global directory_name

load([directory_name '/xgrid']);
load([directory_name '/ygrid']);
[xloc yloc]=meshgrid(xcoor,ycoor);
velo_x1=min(xcoor);
velo_x2=max(xcoor);
velo_y1=min(ycoor);
velo_y2=max(ycoor);

velo_domain=['[' num2str(velo_x1,'%8.2f') ', ' num2str(velo_x2,'%8.2f')...
                ', ' num2str(velo_y1,'%8.2f') ', ' num2str(velo_y2,'%8.2f') ']'];
helpdlg(['Please set the parameters for calculation region and meshsize. ' ...
            'The calculation region [X1, X2, Y1, Y2] must be a subdomain of the velocity field'...
            velo_domain '.'],...
            'Calculation Region')



% --------------------------------------------------------------------
function ComputeFTLEMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to ComputeFTLEMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global nx ny f_start f_end f_int

f_start=str2double(get(handles.starting_frame,'String'));
f_end=str2double(get(handles.end_frame,'String'));
f_int=str2double(get(handles.interval,'String'));
t_length=str2double(get(handles.integration_length,'String'));
tstep=str2double(get(handles.time_interval,'String'));
dt=str2double(get(handles.step_size,'String'));

region_x1=str2double(get(handles.region_x_min,'String'));
region_x2=str2double(get(handles.region_x_max,'String'));
region_y1=str2double(get(handles.region_y_min,'String'));
region_y2=str2double(get(handles.region_y_max,'String'));

nx=str2double(get(handles.meshsize_x,'String'));
ny=str2double(get(handles.meshsize_y,'String'));

forward_calculation=double(get(handles.forward_FTLE,'Value'));
backward_calculation=double(get(handles.backward_FTLE,'Value'));


% f_start
% t_length
% tstep
% dt
% region_x1
% region_x2
% region_y1
% region_y2
% nx
% ny
% forward_calculation
% backward_calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%some variables repeated recalled by subroutines are set as global
%variables to save computation time
global xloc yloc xcoor ycoor sxcoor sycoor umesh vmesh directory_name sigma

load([directory_name '/xgrid']);
load([directory_name '/ygrid']);
[xloc yloc]=meshgrid(xcoor,ycoor);
velo_x1=min(xcoor);
velo_x2=max(xcoor);
velo_y1=min(ycoor);
velo_y2=max(ycoor);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sxcoor=linspace(region_x1,region_x2,nx);               % define the initial region of the fluid to track
sycoor=linspace(region_y1,region_y2,ny);         
[xmesh ymesh]=meshgrid(sxcoor,sycoor);
sigma_x=xmesh';                         
sigma_y=ymesh';                         
sigma=zeros(nx,ny);
sigma0=sigma;

CalcFTLE=zeros(nx,ny)+1;                
LeftDomain=zeros(nx,ny);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Forward FTLE Calculation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if forward_calculation==1
    for iT=f_start:f_int:f_end                        % define the time when the FTLE is calculated

        % calculate FTLE at frame # iT
        t_proc = ['T = ', ...
            num2str( iT )];
        disp( t_proc );

        t_start=iT;                         % starting frame
        t_end=t_start+t_length;
        T_start=(t_start-1)*tstep;          % initial time
        T_end=(t_end-1)*tstep;              % end time
        T_span=T_end-T_start;               % integration time length

        sigma_x=xmesh';                     % initial position
        sigma_y=ymesh';

        CalcFTLE=zeros(nx,ny)+1;
        LeftDomain=zeros(nx,ny);

        for t_integration=0:dt:t_length-dt;
            %integrate the flow trajectory from t0 to t1;
            t0=t_start+t_integration;
            t1=t0+dt;

            iiT=t0;
            %%%%%%%%%%%% read the velocity files
            expression=['load(''' directory_name '/U_T' num2str(iiT) '.mat'');'];
            eval(expression);
            expression=['load(''' directory_name '/V_T' num2str(iiT) '.mat'');'];
            eval(expression);
            expression=['umesh=u' num2str(iiT) ';'];
            eval(expression);
            expression=['vmesh=v' num2str(iiT) ';'];
            eval(expression);
            expression=['clear u' num2str(iiT) ' v' num2str(iiT) ';'];
            eval(expression);

            clear X Y T
            %%%%%% if the point is inside, calculate it's trajectory
            [X,Y,T]= RungeKutta_f(t0,t1,tstep,sigma_x,sigma_y);
            sigma_x_new=X;
            sigma_y_new=Y;
            
            for ix=1:nx
                for iy=1:ny
                    clear X Y T
                    %%%%%% if the point is inside, calculate it's trajectory
                    if LeftDomain(ix,iy)==0;

                        % of the point is convected out of the domian,
                        % calculate the FTLE at this point and adjacent points.
                        if (sigma_x_new(ix,iy)-velo_x1)*(sigma_x_new(ix,iy)-velo_x2)>=0 | ...
                                (sigma_y_new(ix,iy)-velo_y1)*(sigma_y_new(ix,iy)-velo_y2)>=0
                            LeftDomain(ix,iy)=1;

                            if CalcFTLE(ix,iy)==1
                                sigma(ix,iy)=calculate_FTLE(ix,iy,sigma_x,sigma_y,T_span);
                                CalcFTLE(ix,iy)=0;

                                if ix>1 & CalcFTLE(ix-1,iy)==1
                                    sigma(ix-1,iy)=calculate_FTLE(ix-1,iy,sigma_x,sigma_y,T_span);
                                    CalcFTLE(ix-1,iy)=0;
                                end
                                if  ix<nx & CalcFTLE(ix+1,iy)==1
                                    sigma(ix+1,iy)=calculate_FTLE(ix+1,iy,sigma_x,sigma_y,T_span);
                                    CalcFTLE(ix+1,iy)=0;
                                end
                                if  iy>1 & CalcFTLE(ix,iy-1)==1
                                    sigma(ix,iy-1)=calculate_FTLE(ix,iy-1,sigma_x,sigma_y,T_span);
                                    CalcFTLE(ix,iy-1)=0;
                                end
                                if  iy<ny & CalcFTLE(ix,iy+1)==1
                                    sigma(ix,iy+1)=calculate_FTLE(ix,iy+1,sigma_x,sigma_y,T_span);
                                    CalcFTLE(ix,iy+1)=0;
                                end

                            end
                        end
                        %%%%%% if the point is outside, don't calculate it's trajectory. keep its position
                    else
                        sigma_x_new(ix,iy)=sigma_x(ix,iy);
                        sigma_y_new(ix,iy)=sigma_y(ix,iy);
                    end
                end
            end

            % update the position for all the points in the domain
            for ix=1:nx
                for iy=1:ny
                    if LeftDomain(ix,iy)==0
                        sigma_x(ix,iy)=sigma_x_new(ix,iy);
                        sigma_y(ix,iy)=sigma_y_new(ix,iy);
                    end
                end
            end
            
            c_frame = (['Frame = ',num2str(iT)]);
            c_proc = strcat( 'process accomplished :  ', ...
                num2str( 100*t_integration/t_length,' %03.0f' ), '/100' );
            disp( c_proc );


            axes(handles.axes1);
            Z = sigma0';
            [C,h] = contourf(xmesh,ymesh,Z,10);
            axis([region_x1 region_x2 region_y1 region_y2]);
            
            xlabel('X');
            ylabel('Y');
            title('FTLE plot');
            text_output=strvcat(c_frame,c_proc);
            text(region_x2,region_y2,text_output);
            pause(0.2);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end

        % calculate FTLE for all the points inside the domain at the last
        % time step.
        c_proc = strcat( 'process accomplished :  ', ...
            num2str( 100,' %03.0f' ), '/100' );

        for ix=1:nx
            for iy=1:ny
                if CalcFTLE(ix,iy)==1
                    sigma(ix,iy)=calculate_FTLE(ix,iy,sigma_x,sigma_y,T_span);
                    CalcFTLE(ix,iy)=0;
                end
            end
        end
        
        sigma0=sigma;
        expression=['clear FTLE' num2str(iT)];
        eval(expression);
        expression=['global FTLE' num2str(iT)];
        eval(expression);
        expression=['FTLE' num2str(iT) '=sigma;'];
        eval(expression);
        
        
        %     figure
        axes(handles.axes1);
        Z = sigma0';
        [C,h] = contourf(xmesh,ymesh,Z,10);
        axis([region_x1 region_x2 region_y1 region_y2]);
        
        xlabel('X');
        ylabel('Y');
        title('FTLE plot');
        text_output=strvcat(c_frame,c_proc);
        text(region_x2,region_y2,text_output);

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Backward FTLE Calculation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if backward_calculation==1
    for iT=f_start:f_int:f_end                        % define the time when the FTLE is calculated

        % calculate FTLE at frame # iT
        t_proc = ['T = ', ...
            num2str( iT )];
        disp( t_proc );

        t_start=iT;                         % starting frame
        t_end=t_start-t_length;
        T_start=(t_start-1)*tstep;          % initial time
        T_end=(t_end-1)*tstep;              % end time
        T_span=abs(T_end-T_start);               % integration time length

        sigma_x=xmesh';                     % initial position
        sigma_y=ymesh';
        
        CalcFTLE=zeros(nx,ny)+1;
        LeftDomain=zeros(nx,ny);

        for t_integration=0:dt:t_length-dt;
            %integrate the flow trajectory from t0 to t1;
            t0=t_start-t_integration;
            t1=t0-dt;

            iiT=t0;
            %%%%%%%%%%%% read the velocity files
            expression=['load(''' directory_name '/U_T' num2str(iiT) '.mat'');'];
            eval(expression);
            expression=['load(''' directory_name '/V_T' num2str(iiT) '.mat'');'];
            eval(expression);
            expression=['umesh=u' num2str(iiT) ';'];
            eval(expression);
            expression=['vmesh=v' num2str(iiT) ';'];
            eval(expression);
            expression=['clear u' num2str(iiT) ' v' num2str(iiT) ';'];
            eval(expression);

            clear X Y T
            %%%%%% if the point is inside, calculate it's trajectory
            [X,Y,T]= RungeKutta_b(t0,t1,tstep,sigma_x,sigma_y);
            sigma_x_new=X;
            sigma_y_new=Y;

            for ix=1:nx
                for iy=1:ny
                    clear X Y T
                    %%%%%% if the point is inside, calculate it's trajectory
                    if LeftDomain(ix,iy)==0;

                        % of the point is convected out of the domian,
                        % calculate the FTLE at this point and adjacent points.
                        if (sigma_x_new(ix,iy)-velo_x1)*(sigma_x_new(ix,iy)-velo_x2)>=0 | ...
                                (sigma_y_new(ix,iy)-velo_y1)*(sigma_y_new(ix,iy)-velo_y2)>=0
                            LeftDomain(ix,iy)=1;

                            if CalcFTLE(ix,iy)==1
                                sigma(ix,iy)=calculate_FTLE(ix,iy,sigma_x,sigma_y,T_span);
                                CalcFTLE(ix,iy)=0;

                                if ix>1 & CalcFTLE(ix-1,iy)==1
                                    sigma(ix-1,iy)=calculate_FTLE(ix-1,iy,sigma_x,sigma_y,T_span);
                                    CalcFTLE(ix-1,iy)=0;
                                end
                                if  ix<nx & CalcFTLE(ix+1,iy)==1
                                    sigma(ix+1,iy)=calculate_FTLE(ix+1,iy,sigma_x,sigma_y,T_span);
                                    CalcFTLE(ix+1,iy)=0;
                                end
                                if  iy>1 & CalcFTLE(ix,iy-1)==1
                                    sigma(ix,iy-1)=calculate_FTLE(ix,iy-1,sigma_x,sigma_y,T_span);
                                    CalcFTLE(ix,iy-1)=0;
                                end
                                if  iy<ny & CalcFTLE(ix,iy+1)==1
                                    sigma(ix,iy+1)=calculate_FTLE(ix,iy+1,sigma_x,sigma_y,T_span);
                                    CalcFTLE(ix,iy+1)=0;
                                end

                            end
                        end
                        %%%%%% if the point is outside, don't calculate it's trajectory. keep its position
                    else
                        sigma_x_new(ix,iy)=sigma_x(ix,iy);
                        sigma_y_new(ix,iy)=sigma_y(ix,iy);
                    end
                end
            end

            % update the position for all the points in the domain
            for ix=1:nx
                for iy=1:ny
                    if LeftDomain(ix,iy)==0
                        sigma_x(ix,iy)=sigma_x_new(ix,iy);
                        sigma_y(ix,iy)=sigma_y_new(ix,iy);
                    end
                end
            end
            
            c_frame = (['Frame = ',num2str(iT)]);
            c_proc = strcat( 'process accomplished :  ', ...
                num2str( 100*t_integration/t_length,' %03.0f' ), '/100' );
            disp( c_proc );

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             for ix=1:nx
%                 for iy=1:ny
%                     if CalcFTLE(ix,iy)==1
%                         sigma(ix,iy)=calculate_FTLE(ix,iy,sigma_x,sigma_y,T_span);
%                         %                         CalcFTLE(ix,iy)=0;
%                     end
%                 end
%             end



            %     figure
            axes(handles.axes1);
            Z = sigma0';
            [C,h] = contourf(xmesh,ymesh,Z,10);
            axis([region_x1 region_x2 region_y1 region_y2]);
            
            xlabel('X');
            ylabel('Y');
            title('FTLE plot');
            text_output=strvcat(c_frame,c_proc);
            text(region_x2,region_y2,text_output);
            pause(0.2);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end

        % calculate FTLE for all the points inside the domain at the last
        % time step.
        c_proc = strcat( 'process accomplished :  ', ...
            num2str( 100,' %03.0f' ), '/100' );

        for ix=1:nx
            for iy=1:ny
                if CalcFTLE(ix,iy)==1
                    sigma(ix,iy)=calculate_FTLE(ix,iy,sigma_x,sigma_y,T_span);
                    CalcFTLE(ix,iy)=0;
                end
            end
        end
        
        sigma0=sigma;
        
        expression=['clear FTLE' num2str(iT)];
        eval(expression);
        expression=['global FTLE' num2str(iT)];
        eval(expression);
        expression=['FTLE' num2str(iT) '=sigma;'];
        eval(expression);
        
        
        %     figure
        axes(handles.axes1);
        Z = sigma0';
        [C,h] = contourf(xmesh,ymesh,Z,10);
        axis([region_x1 region_x2 region_y1 region_y2]);
        
        xlabel('X');
        ylabel('Y');
        title('FTLE plot');
        text_output=strvcat(c_frame,c_proc);
        text(region_x2,region_y2,text_output);

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Menu Control.4.Output%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------------------
function OutputMenu_Callback(hObject, eventdata, handles)
% hObject    handle to OutputMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function SaveFTLECalculationMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to SaveFTLECalculationMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global sigma
[output_file,output_path] = uiputfile('FTLE.txt','Save file name');
output_FTLE(sigma,output_file,output_path)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Menu Control.5.Help%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------------------
function HelpMenu_Callback(hObject, eventdata, handles)
% hObject    handle to HelpMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
% --------------------------------------------------------------------
function ReadMeMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to ReadMeMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

open('ReadMe.pdf');
% --------------------------------------------------------------------
function AboutMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to AboutMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
message{1}=['FTLE Matlab Kit V1.0'];
message{2}=['California Institute of Technology']; 
msgbox(message,'About')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function sort_char_array%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sortedArr,finalIndex] = sort_char_array(orgArr)

% SORTN - Sort cell array or char array in textual and numerical mode.
%               Primary sort by first letter, secondary by second letter and so on.
%               Numbers are considered an numerical values.
%               
%               Sorting Rules:
%                           'a' comes before 'b' and so on.
%                           Numbers come before letters
%                           Numbers are compared according to their value
%                           Short string comes before a long one (when first chars are identical) 
%                           Ignore characters that are not letters or
%                           numbers.
%                           SORTN is NOT case sensitive.
%
%               The main interest in this function is that it can sort
%               numbers with text, so string like 'a60' comes before 'a100'
%
%[SORTED_CELL_ARRAY , INDEXES] = SORTN(CELL_ARRAY) - Sorts the cells array
%                or char array and returns sorted cell array and indexes
%                vector of new arranged cell array.
%
%EXAMPLE:
%               cell_arr = {'abcde54f';
%                                       'aabde54f'}
%
%              [sorted,indexes]=sortn(cell_arr)
%
%               sorted = 
%                                   'aabde54f'
%                                   'abcde54f'
%               indexes = 
%                                       2
%                                       1                          
%
%
%
%               cell_arr = {'abcde154f';
%                                       'abcde54f'}
%
%              [sorted,indexes]=sortn(cell_arr)
%
%               sorted = 
%
%                                   'abcde54f'
%                                   'abcde154f'
%               indexes =
%
%                                   2
%                                   1    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if ischar(orgArr)   %function works with cell array, so if char input, convert to cell array
    tempList = lower(orgArr);
    orgArr = cellstr(orgArr);
else
    tempList = char(lower(orgArr));
end

finalIndex = (1:length(orgArr))';   %create the index vector to be changed
g = {};
for i=1:size(tempList)  %go over all strings and break each one to single letters and numbers
    tempCellArr = breakString(tempList(i,:));   %break one string from the cell array
    g = AssignCellArray(tempCellArr, g);    %add the breaked string to cell matrix
end
%[g,indexes] = SortByColumn2(g,5);
for i=size(g,2):-1:1    %go over cell matrix of breaked strings and sort is from last column to the first

    [g,indexes] = SortByColumn(g,i);    %sort the i column
    
    orgArr = orgArr(indexes);   %update cell array order by the indexes
    finalIndex = finalIndex(indexes);   %update final index order
    
end

sortedArr =orgArr ;


function [sortedArr,indexes] = SortByColumn(arr,column);    
%this function sort each column of the breaked cell-matrix
%sorting is done by converting each textual number to its real value
%and later assign every letter its ascii value plus the maximum of all
%numbers. 
%the sorting is pure numeric (Matlab sort) so numbers are always before
%letters.

tempColumn = zeros(size(arr,1),1);
for i=1:length(tempColumn)
    if ~isempty(str2num(char(arr(i,column)))); %finds all numbers including non-real numbers
        tempColumn(i) = str2num(char(arr(i,column)))+1; %assign numbers to the vactor
    end
end
maxNum = max(tempColumn);
u2 = isletter(char(arr(:,column))); %find all letters
tempColumn(u2) = char(arr(u2,column)) + maxNum; %assign numerical representation of letters that is bigger then the numbers

[tempSort,indexes] = sort(tempColumn);
sortedArr = arr(indexes,:);



function newArr = AssignCellArray(assignmentVector, assignTo)
%this function helps adding cell-arrays to the breaked strings cell-matrix


assignToSize = size(assignTo,2);
assignmentVectorSize = length(assignmentVector);
if assignToSize==0
    newArr = assignmentVector;
    
elseif (assignToSize==assignmentVectorSize)
    newArr = [assignTo;assignmentVector];
elseif (assignToSize>assignmentVectorSize)
    assignmentVector((assignmentVectorSize+1):assignToSize)={' '};
    newArr = [assignTo;assignmentVector];
else
    assignTo(:,(assignToSize+1):assignmentVectorSize)={' '};
    newArr = [assignTo;assignmentVector];
end


function groups = breakString(stringToBreak)
%this function breaks a string to its single letters and numbers,
%producing cell array for a single string.

if ~ischar(stringToBreak)   %convert to char array
    tempString = char(stringToBreak);
else
    tempString=stringToBreak;
end
counter = 1;
i=1;
groups = {};
while i<=length(tempString) %go over all chars in the string
    if  ~isletter(tempString(i)) & isempty(str2num(tempString(i)))   %tempString(i)==' '   %Skip spaces
        i=i+1;
        continue;
    end
    if isempty(str2num(tempString(i)))  %in case of a letter - just add it to the array
        groups(counter)={tempString(i)};
        counter = counter + 1;
        i = i + 1;
    else
        groups(counter) = {tempString(i)};
        i = i+1;
        while (i<=length(tempString) & ~isempty(str2num(tempString(i))))    %in case of number,
            %add it to the array and search for following numbers
            groups(counter)={[char(groups(counter)) tempString(i)]};
            i = i + 1;
        end
        counter = counter + 1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function velocity_plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=velocity_plot(pathname,iT1,iT2,plotwindow)

load([pathname '/xgrid']);
load([pathname '/ygrid']);
[xloc yloc]=meshgrid(xcoor,ycoor);

% figure
axes(plotwindow);
for i=iT1:iT2
    expression=['load(''' pathname '/U_T' num2str(i) '.mat'');'];
    eval(expression);
    expression=['load(''' pathname '/V_T' num2str(i) '.mat'');'];
    eval(expression);
    expression=['u=u' num2str(i) ';'];
    eval(expression);
    expression=['v=v' num2str(i) ';'];
    eval(expression);
    expression=['clear u' num2str(i) ' v' num2str(i) ';'];
    eval(expression);
    clear load_UX load_UY;
    
    quiver(xloc,yloc,u',v',1);
    axis([min(xcoor) max(xcoor) min(ycoor) max(ycoor)]);
    %     axis([-12 2 -3 3]);
    xlabel('X');
    ylabel('Y');
    title('velocity plot');
    text1=['slide # :' num2str(i)];
    text(max(xcoor)-1,max(ycoor)+.5,text1);
    pause(0.2);
    out=i;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function RungeKutta_f%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X,Y,T]= RungeKutta_f(i_start,i_end,h,xstart,ystart)

% global umesh vmesh

n=1;
t=i_start;
tend=i_end;

delta_i=i_end-i_start;

[ux,uy] = velocity(xstart,ystart,i_start);
x = xstart + delta_i*h*ux;
y = ystart + delta_i*h*uy;
    
X=x;
Y=y;
T=t+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function RungeKutta_b%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Y,T]= RungeKutta_b(i_start,i_end,h,xstart,ystart)

% global umesh vmesh

n=1;
t=i_start;
tend=i_end;

delta_i=abs(i_end-i_start);

[ux,uy] = velocity(xstart,ystart,i_start);
x = xstart - delta_i*h*ux;
y = ystart - delta_i*h*uy;

X=x;
Y=y;
T=t-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function velocity%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ux,uy] = velocity(x,y,iiT)

global xloc yloc umesh vmesh
    
ux = interp2(xloc,yloc,umesh',x,y,'linear');
uy = interp2(xloc,yloc,vmesh',x,y,'linear');

index=find(isnan(ux)==1);
ux(index)=0;

index=find(isnan(uy)==1);
uy(index)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function calculate_FTLE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out=calculate_FTLE(ix,iy,flowmap_x,flowmap_y,t_span)

global sxcoor sycoor

[nx ny]=size(flowmap_x);

if (ix-1)*(ix-nx)<0 & (iy-1)*(iy-ny)<0
    A11=(flowmap_x(ix+1,iy)-flowmap_x(ix-1,iy))/(sxcoor(ix+1)-sxcoor(ix-1));
    A12=(flowmap_x(ix,iy+1)-flowmap_x(ix,iy-1))/(sycoor(iy+1)-sycoor(iy-1));
    A21=(flowmap_y(ix+1,iy)-flowmap_y(ix-1,iy))/(sxcoor(ix+1)-sxcoor(ix-1));
    A22=(flowmap_y(ix,iy+1)-flowmap_y(ix,iy-1))/(sycoor(iy+1)-sycoor(iy-1));
    A=[A11 A12;A21 A22];
    B=A'*A;
    delta=max(eig(B));
    out=log(delta)/(2*t_span);
else
    out=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function output%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=output_FTLE(sigma,file,path)

global sxcoor sycoor f_start f_end f_int

% f_start=str2double(get(handles.starting_frame,'String'));
% f_end=str2double(get(handles.end_frame,'String'));
% f_int=str2double(get(handles.interval,'String'));

%     expression=['save(''sigma_f'',''sigma'');'];
%     eval(expression);
    
[xmesh ymesh]=meshgrid(sxcoor,sycoor);
xloc1=reshape(xmesh(:,:),[],1);
yloc1=reshape(ymesh(:,:),[],1);

filename=[path file];
fid = fopen(filename,'wt');
fprintf(fid,'TITLE = FTLE FIELD\n');
fprintf(fid,'VARIABLES = X, Y, FTLE\n');
for iT=f_start:f_int:f_end
    expression=['global FTLE' num2str(iT) ';'];
    eval(expression);
    expression=['sigma0=FTLE' num2str(iT) ';'];
    eval(expression);

    [nx,ny]=size(sigma0);

    % FTLE FIELD
    Z = sigma0';

    Z1=reshape(Z,[],1);
    dummy=0*Z1;
    data=[xloc1 yloc1 Z1];

    fprintf(fid,['ZONE T=''A' num2str(iT) ''' I=' num2str(ny) ', J=' num2str(nx) ', F=POINT\n']);
    fprintf(fid,'%12.8f %12.8f %12.8f \n',data');
end
fclose(fid);















% --------------------------------------------------------------------
function ParticleFTLECalculationMenu_Callback(hObject, eventdata, handles)
% hObject    handle to ParticleFTLECalculationMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ParametersMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to ParametersMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mass_ratio inertia

prompt = {'Mass ratio parameter R:','Inertia parameter A:'};
dlg_title = 'Enter values for Maxey-Riley equation parameters ';
num_lines = 1;
def = {'2/3','2'};
answer = inputdlg(prompt,dlg_title,num_lines,def);

mass_ratio=str2num(answer{1});
inertia=str2num(answer{2});

% --------------------------------------------------------------------
function ComputeParticleFTLEMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to ComputeParticleFTLEMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global nx ny f_start f_end f_int mass_ratio inertia

f_start=str2double(get(handles.starting_frame,'String'));
f_end=str2double(get(handles.end_frame,'String'));
f_int=str2double(get(handles.interval,'String'));
t_length=str2double(get(handles.integration_length,'String'));
tstep=str2double(get(handles.time_interval,'String'));
dt=str2double(get(handles.step_size,'String'));

region_x1=str2double(get(handles.region_x_min,'String'));
region_x2=str2double(get(handles.region_x_max,'String'));
region_y1=str2double(get(handles.region_y_min,'String'));
region_y2=str2double(get(handles.region_y_max,'String'));

nx=str2double(get(handles.meshsize_x,'String'));
ny=str2double(get(handles.meshsize_y,'String'));

% f_start
% t_length
% tstep
% dt
% region_x1
% region_x2
% region_y1
% region_y2
% nx
% ny
% forward_calculation
% backward_calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%some variables repeated recalled by subroutines are set as global
%variables to save computation time
global xloc yloc xcoor ycoor sxcoor sycoor umesh vmesh directory_name sigma

load([directory_name '/xgrid']);
load([directory_name '/ygrid']);
[xloc yloc]=meshgrid(xcoor,ycoor);
velo_x1=min(xcoor);
velo_x2=max(xcoor);
velo_y1=min(ycoor);
velo_y2=max(ycoor);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sxcoor=linspace(region_x1,region_x2,nx);               % define the initial region of the fluid to track
sycoor=linspace(region_y1,region_y2,ny);         
[xmesh ymesh]=meshgrid(sxcoor,sycoor);
sigma_x=xmesh';                         
sigma_y=ymesh';                         
sigma=zeros(nx,ny);
sigma0=sigma;

R=mass_ratio;
A=inertia;
CalcFTLE=zeros(nx,ny)+1;                
LeftDomain=zeros(nx,ny);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Forward FTLE Calculation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iT=f_start:f_int:f_end                        % define the time when the FTLE is calculated

    % calculate FTLE at frame # iT
    t_proc = ['T = ', ...
        num2str( iT )];
    disp( t_proc );

    t_start=iT;                         % starting frame
    t_end=t_start+t_length;
    T_start=(t_start-1)*tstep;          % initial time
    T_end=(t_end-1)*tstep;              % end time
    T_span=T_end-T_start;               % integration time length
    delta_t=tstep*dt;
    
%     delta_x=(xcoor(2)-xcoor(1))/5;
%     delta_y=(ycoor(2)-ycoor(1))/5;
    delta_x=1;
    delta_y=1;
    
    sigma_x=xmesh';                     % initial position
    sigma_y=ymesh';

    CalcFTLE=zeros(nx,ny)+1;
    LeftDomain=zeros(nx,ny);
    
    %initial fluid velocity
    expression=['load(''' directory_name '/U_T' num2str(t_start) '.mat'');'];
    eval(expression);
    expression=['load(''' directory_name '/V_T' num2str(t_start) '.mat'');'];
    eval(expression);
    expression=['umesh=u' num2str(t_start) ';'];
    eval(expression);
    expression=['vmesh=v' num2str(t_start) ';'];
    eval(expression);
    expression=['clear u' num2str(t_start) ' v' num2str(t_start) ';'];
    eval(expression);

    %initial particle velocity
    [vx,vy]=velocity(sigma_x,sigma_y,t_start);
    
    for t_integration=0:dt:t_length-dt;
        %integrate the flow trajectory from t0 to t1;
        t0=t_start+t_integration;
%         t1=t0+dt;

        iiT=t0;
        %%%%%%%%%%%% read the velocity files
        expression=['load(''' directory_name '/U_T' num2str(iiT) '.mat'');'];
        eval(expression);
        expression=['load(''' directory_name '/V_T' num2str(iiT) '.mat'');'];
        eval(expression);
        expression=['umesh=u' num2str(iiT) ';'];
        eval(expression);
        expression=['vmesh=v' num2str(iiT) ';'];
        eval(expression);
        expression=['clear u' num2str(iiT) ' v' num2str(iiT) ';'];
        eval(expression);
        
        
        %fluid velocity
        [ux,uy]=velocity(sigma_x,sigma_y,iiT);
        
        %fluid velocity at neighbor points
        [ux_e,uy_e]=velocity(sigma_x+delta_x,sigma_y,iiT);
        [ux_w,uy_w]=velocity(sigma_x-delta_x,sigma_y,iiT);
        [ux_n,uy_n]=velocity(sigma_x,sigma_y+delta_y,iiT);
        [ux_s,uy_s]=velocity(sigma_x,sigma_y-delta_y,iiT);
        
        %fluid velocity at previous frame
        expression=['load(''' directory_name '/U_T' num2str(iiT-1) '.mat'');'];
        eval(expression);
        expression=['load(''' directory_name '/V_T' num2str(iiT-1) '.mat'');'];
        eval(expression);
        expression=['umesh=u' num2str(iiT-1) ';'];
        eval(expression);
        expression=['vmesh=v' num2str(iiT-1) ';'];
        eval(expression);
        expression=['clear u' num2str(iiT-1) ' v' num2str(iiT-1) ';'];
        eval(expression);
        
        [ux0,uy0]=velocity(sigma_x,sigma_y,iiT-1);
        
        %fluid velocity at next frame
        expression=['load(''' directory_name '/U_T' num2str(iiT+1) '.mat'');'];
        eval(expression);
        expression=['load(''' directory_name '/V_T' num2str(iiT+1) '.mat'');'];
        eval(expression);
        expression=['umesh=u' num2str(iiT+1) ';'];
        eval(expression);
        expression=['vmesh=v' num2str(iiT+1) ';'];
        eval(expression);
        expression=['clear u' num2str(iiT+1) ' v' num2str(iiT+1) ';'];
        eval(expression);

        [ux2,uy2]=velocity(sigma_x,sigma_y,iiT+1);
        
        %derivatives of fluid velocity
        duxdx=(ux_e-ux_w)/(2*delta_x);
        duydx=(uy_e-uy_w)/(2*delta_x);
        duxdy=(ux_n-ux_s)/(2*delta_y);
        duydy=(uy_n-uy_s)/(2*delta_y);
        duxdt=(ux2-ux0)/(2*delta_t);
        duydt=(uy2-uy0)/(2*delta_t);
        
        %Lagrangian acceleration of fluid particles
        DuxDt=duxdt+ux.*duxdx+uy.*duxdy;
        DuyDt=duydt+ux.*duydx+uy.*duydy;
        
        %Lagrangian acceleration of inertial particles
        DvxDt=3*R*DuxDt/2-A*(vx-ux);
        DvyDt=3*R*DuyDt/2-A*(vy-uy);

        %next position of particle
        sigma_x_new=sigma_x+vx*delta_t;
        sigma_y_new=sigma_y+vy*delta_t;

        %velocity of particle at next position
        vx=vx+DvxDt*delta_t;
        vy=vy+DvyDt*delta_t;

        %calculate FTLE
        for ix=1:nx
            for iy=1:ny
                clear X Y T
                %%%%%% if the point is inside, calculate it's trajectory
                if LeftDomain(ix,iy)==0;
                    % of the point is convected out of the domian,
                    % calculate the FTLE at this point and adjacent points.
                    if (sigma_x_new(ix,iy)-velo_x1)*(sigma_x_new(ix,iy)-velo_x2)>=0 | ...
                            (sigma_y_new(ix,iy)-velo_y1)*(sigma_y_new(ix,iy)-velo_y2)>=0
                        LeftDomain(ix,iy)=1;

                        if CalcFTLE(ix,iy)==1
                            sigma(ix,iy)=calculate_FTLE(ix,iy,sigma_x,sigma_y,T_span);
                            CalcFTLE(ix,iy)=0;

                            if ix>1 & CalcFTLE(ix-1,iy)==1
                                sigma(ix-1,iy)=calculate_FTLE(ix-1,iy,sigma_x,sigma_y,T_span);
                                CalcFTLE(ix-1,iy)=0;
                            end
                            if  ix<nx & CalcFTLE(ix+1,iy)==1
                                sigma(ix+1,iy)=calculate_FTLE(ix+1,iy,sigma_x,sigma_y,T_span);
                                CalcFTLE(ix+1,iy)=0;
                            end
                            if  iy>1 & CalcFTLE(ix,iy-1)==1
                                sigma(ix,iy-1)=calculate_FTLE(ix,iy-1,sigma_x,sigma_y,T_span);
                                CalcFTLE(ix,iy-1)=0;
                            end
                            if  iy<ny & CalcFTLE(ix,iy+1)==1
                                sigma(ix,iy+1)=calculate_FTLE(ix,iy+1,sigma_x,sigma_y,T_span);
                                CalcFTLE(ix,iy+1)=0;
                            end

                        end
                    end
                    %%%%%% if the point is outside, don't calculate it's trajectory. keep its position
                else
                    sigma_x_new(ix,iy)=sigma_x(ix,iy);
                    sigma_y_new(ix,iy)=sigma_y(ix,iy);
                end
            end
        end

        % update the position for all the points in the domain
        for ix=1:nx
            for iy=1:ny
                if LeftDomain(ix,iy)==0
                    sigma_x(ix,iy)=sigma_x_new(ix,iy);
                    sigma_y(ix,iy)=sigma_y_new(ix,iy);
                end
            end
        end

        c_frame = (['Frame = ',num2str(iT)]);
        c_proc = strcat( 'process accomplished :  ', ...
            num2str( 100*t_integration/t_length,' %03.0f' ), '/100' );
        disp( c_proc );

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %             for ix=1:nx
        %                 for iy=1:ny
        %                     if CalcFTLE(ix,iy)==1
        %                         sigma(ix,iy)=calculate_FTLE(ix,iy,sigma_x,sigma_y,T_span);
        %                         %                         CalcFTLE(ix,iy)=0;
        %                     end
        %                 end
        %             end
        %
        %
        %
        %             %     figure
        axes(handles.axes1);
        Z = sigma0';
        [C,h] = contourf(xmesh,ymesh,Z,10);
        axis([region_x1 region_x2 region_y1 region_y2]);

        xlabel('X');
        ylabel('Y');
        title('FTLE plot');
        text_output=strvcat(c_frame,c_proc);
        text(region_x2,region_y2,text_output);
        pause(0.2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    % calculate FTLE for all the points inside the domain at the last
    % time step.
    c_proc = strcat( 'process accomplished :  ', ...
        num2str( 100,' %03.0f' ), '/100' );

    for ix=1:nx
        for iy=1:ny
            if CalcFTLE(ix,iy)==1
                sigma(ix,iy)=calculate_FTLE(ix,iy,sigma_x,sigma_y,T_span);
                CalcFTLE(ix,iy)=0;
            end
        end
    end

    sigma0=sigma;
    expression=['clear FTLE' num2str(iT)];
    eval(expression);
    expression=['global FTLE' num2str(iT)];
    eval(expression);
    expression=['FTLE' num2str(iT) '=sigma;'];
    eval(expression);


    %     figure
    axes(handles.axes1);
    Z = sigma0';
    [C,h] = contourf(xmesh,ymesh,Z,10);
    axis([region_x1 region_x2 region_y1 region_y2]);

    xlabel('X');
    ylabel('Y');
    title('FTLE plot');
    text_output=strvcat(c_frame,c_proc);
    text(region_x2,region_y2,text_output);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

