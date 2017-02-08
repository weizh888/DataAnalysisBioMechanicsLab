function varargout = step1_correct_dist_GUI(varargin)
% STEP1_CORRECT_DIST_GUI MATLAB code for step1_correct_dist_GUI.fig
%      STEP1_CORRECT_DIST_GUI, by itself, creates a new STEP1_CORRECT_DIST_GUI or raises the existing
%      singleton*.
%
%      H = STEP1_CORRECT_DIST_GUI returns the handle to a new STEP1_CORRECT_DIST_GUI or the handle to
%      the existing singleton*.
%
%      STEP1_CORRECT_DIST_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STEP1_CORRECT_DIST_GUI.M with the given input arguments.
%
%      STEP1_CORRECT_DIST_GUI('Property','Value',...) creates a new STEP1_CORRECT_DIST_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before step1_correct_dist_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to step1_correct_dist_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help step1_correct_dist_GUI

% Last Modified by GUIDE v2.5 24-Oct-2016 19:33:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @step1_correct_dist_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @step1_correct_dist_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
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


% --- Executes just before step1_correct_dist_GUI is made visible.
function step1_correct_dist_GUI_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to step1_correct_dist_GUI (see VARARGIN)

% Choose default command line output for step1_correct_dist_GUI
handles.output = hObject;

SpeedValue = evalin('base', 'SpeedValue');
Time = evalin('base', 'Time');
A_dist_Y = evalin('base', 'A_dist_Y');
B_dist_Y = evalin('base', 'B_dist_Y');
AB_dist_Y = (A_dist_Y + B_dist_Y)./2.0;

plot(Time, AB_dist_Y, '-r', Time, A_dist_Y, '-b', Time, B_dist_Y, '-g')
legend('(Trap A + B)/2','Trap A','Trap B')
xlabel('Time (s)'); ylabel('Distance (nm)'); title('Distance vs Time')
datacursormode on

set(handles.speed_value, 'string', num2str(SpeedValue))
handles.SpeedValue = SpeedValue;
handles.DistCorrectFactor = [];
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes step1_correct_dist_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = step1_correct_dist_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in point1.
function point1_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to point1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% data_retract = handles.data_retract;
dcm_obj = datacursormode(gcf);
info = getCursorInfo(dcm_obj);
handles.p1 = info.Position;

guidata(hObject, handles);


% --- Executes on button press in point2.
function point2_Callback(hObject, eventdata, handles)
% hObject    handle to point2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dcm_obj = datacursormode(gcf);
info = getCursorInfo(dcm_obj);
handles.p2 = info.Position;

guidata(hObject, handles);

% --- Executes on button press in linear_fit.
function linear_fit_Callback(hObject, eventdata, handles)
% hObject    handle to linear_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Time = evalin('base', 'Time');
A_dist_Y = evalin('base', 'A_dist_Y');
B_dist_Y = evalin('base', 'B_dist_Y');
AB_dist_Y = (A_dist_Y + B_dist_Y)./2.0;

switch get(get(handles.trap_selection,'SelectedObject'),'Tag')
    case 'TrapA_B',  Y_dist = AB_dist_Y;
    case 'TrapA',  Y_dist = A_dist_Y;
    case 'TrapB',  Y_dist = B_dist_Y;
end

index1 = find(Time==handles.p1(1,1));
index2 = find(Time==handles.p2(1,1));
data_x = Time(index1:index2,:);
data_y = Y_dist(index1:index2,:);
       
p = polyfit(data_x,data_y,1);
% xx =
% min(data_x)-0.2*peak2peak(data_x):0.5:max(data_x)+0.2*peak2peak(data_x);
% %peak2peak is defined in signal Processing Toolbox
exPercent = 0.2;
spacing = 0.5;
xx = min(data_x)-exPercent*(max(data_x)-min(data_x)):spacing:max(data_x)+exPercent*(max(data_x)-min(data_x));
yy = p(1)*xx + p(2);
axes(handles.axes1)
hold on
plot(xx,yy,'-m')
hold off

set(handles.slope,'string',num2str(p(1)))
set(handles.dist_factor,'string',num2str(abs(handles.SpeedValue/p(1))))

handles.DistCorrectFactor = [handles.DistCorrectFactor; abs(handles.SpeedValue/p(1))];
set(handles.mean_dist_factor,'string',num2str(mean(handles.DistCorrectFactor)))

guidata(hObject, handles);


function slope_Callback(hObject, eventdata, handles)
% hObject    handle to slope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function slope_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dist_factor_Callback(hObject, eventdata, handles)
% hObject    handle to dist_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function dist_factor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dist_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function speed_value_Callback(hObject, eventdata, handles)
% hObject    handle to speed_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.SpeedValue = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function speed_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to speed_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in done_gui.
function done_gui_Callback(hObject, eventdata, handles)
% hObject    handle to done_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global filename file_path ExpDate
Time = evalin('base', 'Time');
A_dist_Y = evalin('base', 'A_dist_Y');
B_dist_Y = evalin('base', 'B_dist_Y');
AB_dist_Y = (A_dist_Y + B_dist_Y)./2.0;

switch get(get(handles.trap_selection,'SelectedObject'),'Tag')
    case 'TrapA_B',  Y_dist = AB_dist_Y;
    case 'TrapA',  Y_dist = A_dist_Y;
    case 'TrapB',  Y_dist = B_dist_Y;
end

Y_force = evalin('base', 'Y_force');
filename = evalin('base', 'filename');
file_path = ['data/' ExpDate '/' filename '/']; % set file path
new_file = [file_path 'All_Time_Dist_Force.txt']; % [time distance force]
dist_factor_file = ['data/' ExpDate '/Distance_Cor_Factor.txt']; % [dataset factor]

DistCorFactor = mean(handles.DistCorrectFactor);
if length(handles.DistCorrectFactor)==1
    DistCorFactorSTD = 0;
else
    DistCorFactorSTD = std(handles.DistCorrectFactor);
end
    
fid0 = fopen(dist_factor_file,'a+');
fprintf(fid0,'%s %.4f %.4f\r\n',filename,DistCorFactor,DistCorFactorSTD);
fclose(fid0);
% correct distance
Y_dist = Y_dist.*DistCorFactor; % unit of nm
%% Save data
TimeDistForce = horzcat(Time, Y_dist, -Y_force);

% TimeDistForce = [Time'; Y_dist'; -Y_force'];
fid1 = fopen(new_file,'w');
fprintf(fid1,'%8.3f %9.3f %6.2f\r\n',TimeDistForce');
fclose(fid1);

assignin('base', 'DistCorFactor', DistCorFactor)
assignin('base', 'TimeDistForce', TimeDistForce')
close(step1_correct_dist_GUI)


function mean_dist_factor_Callback(hObject, eventdata, handles)
% hObject    handle to mean_dist_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.DistCorrectFactor = str2double(get(handles.mean_dist_factor,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function mean_dist_factor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mean_dist_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
