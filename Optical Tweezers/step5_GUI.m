function varargout = step5_GUI(varargin)
% STEP5_GUI MATLAB code for step5_GUI.fig
%      STEP5_GUI, by itself, creates a new STEP5_GUI or raises the existing
%      singleton*.
%
%      H = STEP5_GUI returns the handle to a new STEP5_GUI or the handle to
%      the existing singleton*.
%
%      STEP5_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STEP5_GUI.M with the given input arguments.
%
%      STEP5_GUI('Property','Value',...) creates a new STEP5_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before step5_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to step5_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help step5_GUI

% Last Modified by GUIDE v2.5 13-Dec-2016 16:43:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @step5_GUI_OpeningFcn, ...
    'gui_OutputFcn',  @step5_GUI_OutputFcn, ...
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


% --- Executes just before step5_GUI is made visible.
function step5_GUI_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to step5_GUI (see VARARGIN)

% Choose default command line output for step5_GUI
handles.output = hObject;

% set( gcf, 'menubar', 'figure' )
clc;
set(handles.axes1,'xtick',[],'ytick',[]);
% load guidata
global kB T filename file_path Index_curves SpeedValue ExpDate
fid = fopen('InfoForStep5.txt','r');
M = textscan(fid,'%s','Delimiter',' ');
ExpDate = M{1}{2};
filename = M{1}{4};
SpeedValue = str2double(M{1}{8});
Index_curves = cell2mat(cellfun(@str2num,M{:}(11:end),'un',0).');
kB = 1.3806503*1e-23; % Boltzmann constant
T = 273.15+24; % absolute temperature
handles.ForceCorFactor = 1;
handles.TargetForce = 67; % DNA OS force, unit: pN
if SpeedValue >= 150
    handles.smoothSpan = 3; % number of points for smoothing
else
    handles.smoothSpan = 7; % number of points for smoothing
end

handles.unfolding_num = 0;
handles.refolding_num = 0;
handles.rupture_num = 0;

handles.shift_true = 1;
handles.shiftDist = 50; % unit: nm
handles.XSpan = 2500;
% SpeedValue = 200;
set(handles.file_name, 'String', filename);
set(handles.speed_value, 'String', SpeedValue);
set(handles.target_force, 'String', handles.TargetForce);
set(handles.smooth_span, 'String', handles.smoothSpan);
set(handles.shift_dist, 'String', handles.shiftDist);

% load list box
file_path = ['data/' ExpDate '/' filename '/'];
source_files = dir(fullfile(file_path, '*.txt'));
single_cycle_file = [ file_path 'Dist_Force_Time ']; % one approaching/retracting process
single_cycle_smooth_file = [ file_path 'DFT ']; % one approaching/retracting process

% check and remove these files not exist in selected curves
for i = 1:length(Index_curves)
    file1 = [single_cycle_file num2str(Index_curves(i)) 'R.txt'];
    file2 = [single_cycle_smooth_file num2str(Index_curves(i)) 'R.txt'];
    if (exist(file1, 'file')~=2) && (exist(file2, 'file')~=2)
        Index_curves(i) = NaN;
    end
end
Index_curves = Index_curves(~isnan(Index_curves));

handles.ForceIndex = 1;
ForceNum = Index_curves(handles.ForceIndex);

str1 = [single_cycle_file num2str(ForceNum) 'A.txt'];
str2 = [single_cycle_file num2str(ForceNum) 'R.txt'];

% check if the data is smoothed or not
flag = 0; % check whether to correct force
if exist(str2, 'file') ~= 2
    str1 = [single_cycle_smooth_file num2str(ForceNum) 'A.txt'];
    str2 = [single_cycle_smooth_file num2str(ForceNum) 'R.txt'];
    flag = 1;
end

handles.data_retract_raw = [];
handles.data_approach_raw = [];
handles.data_retract_cor = [];
handles.data_approach_cor = [];
handles.data_retract = [];
handles.data_approach = [];

x_indent = 0;
y_indent = 0;

handles.data_retract_raw = dlmread(str2);
extension_R_cor = handles.data_retract_raw(:,1) - x_indent; % unit: nm
if flag==1
    force_R_cor = (handles.data_retract_raw(:,2) - y_indent);
else
    force_R_cor = (handles.data_retract_raw(:,2) - y_indent)*handles.ForceCorFactor; % unit: pN
end
time_R_raw = handles.data_retract_raw(:,3); % unit:s for calculating loading rate dF/dt
handles.data_retract_cor = [extension_R_cor force_R_cor time_R_raw];
handles.data_retract = handles.data_retract_cor; % for smoothing data

axes(handles.axes1); cla reset;
% plot retract curve
plot(handles.data_retract_cor(:,1),handles.data_retract_cor(:,2),'-r')

% check if approach file exists: if it does, then plot it
handles.ApproachExist = 0;
if exist(str1, 'file') == 2 %
    handles.ApproachExist = 1;
    handles.data_approach_raw = dlmread(str1);
    extension_A_cor = handles.data_approach_raw(:,1) - x_indent; % unit: nm
    if flag==1
        force_A_cor = (handles.data_approach_raw(:,2) - y_indent); % unit: pN
    else
        force_A_cor = (handles.data_approach_raw(:,2) - y_indent)*handles.ForceCorFactor; % unit: pN
    end
    time_A_raw = handles.data_approach_raw(:,3); % unit:s for calculating loading rate dF/dt
    handles.data_approach_cor = [extension_A_cor force_A_cor time_A_raw];
    handles.data_approach = handles.data_approach_cor;
    
    hold on
    if handles.shift_true==1
        handles.h2 = plot(handles.data_approach_cor(:,1)+handles.shiftDist, handles.data_approach_cor(:,2),'-b');
    else
        handles.h2 = plot(handles.data_approach_cor(:,1), handles.data_approach_cor(:,2),'-b');
    end
    hold off
end
xlabel('Extension (nm)');ylabel('Force (pN)');
title(['Force vs Extension (# ' num2str(Index_curves(handles.ForceIndex)) ')'])
% grid on
% hdt = datacursormode;
% set(hdt,'DisplayStyle','window');
ylim([-5 75]); % set ylim

xLimits = get(gca,'XLim'); yLimits = get(gca,'YLim');
xLimits(2) = xLimits(1)+handles.XSpan;

set(handles.x_min,'string',num2str(xLimits(1),'%5d'));
set(handles.x_max,'string',num2str(xLimits(2),'%5d'));
set(handles.y_min,'string',num2str(yLimits(1),'%5d'));
set(handles.y_max,'string',num2str(yLimits(2),'%5d'));
set(handles.x_span,'string',num2str(xLimits(2)-xLimits(1),'%5d'));
set(handles.y_span,'string',num2str(yLimits(2)-yLimits(1),'%5d'));

handles.XSpan = xLimits(2)-xLimits(1);
handles.YSpan = yLimits(2)-yLimits(1);
handles.axisRange = [xLimits yLimits];
handles.axisRangeZoomed = [xLimits yLimits];
handles.source_files = source_files;
handles.Index_curves = Index_curves;
handles.SpeedValue = SpeedValue;
handles.ExpDate = ExpDate;
% update_range handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = step5_GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function speed_value_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to speed_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SpeedValue = str2double(get(handles.speed_value,'String'));
set(handles.speed_value,'String',num2str(SpeedValue));

handles.SpeedValue = SpeedValue;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function speed_value_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to speed_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in next_file.
function next_file_Callback(hObject, eventdata, handles)
% hObject    handle to next_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global file_path Index_curves
Index_curves = handles.Index_curves;
if handles.ForceIndex<length(Index_curves) && length(Index_curves)>1
    % read data file(s)
    handles.ForceIndex = handles.ForceIndex + 1;
else
    handles.ForceIndex = handles.ForceIndex;
end
ForceNum = Index_curves(handles.ForceIndex);

single_cycle_file = [file_path 'Dist_Force_Time '];
single_cycle_smooth_file = [file_path 'DFT '];

str1 = [single_cycle_file num2str(ForceNum) 'A.txt'];
str2 = [single_cycle_file num2str(ForceNum) 'R.txt'];

% check if the data is smoothed or not
flag = 0; % check whether to correct force
if exist(str2, 'file') ~= 2
    str1 = [single_cycle_smooth_file num2str(ForceNum) 'A.txt'];
    str2 = [single_cycle_smooth_file num2str(ForceNum) 'R.txt'];
    flag = 1;
end

handles.data_retract_raw = [];
handles.data_approach_raw = [];
handles.data_retract_cor = [];
handles.data_approach_cor = [];
handles.data_retract = [];
handles.data_approach = [];

x_indent = 0;
y_indent = 0;

handles.data_retract_raw = dlmread(str2);
extension_R_cor = handles.data_retract_raw(:,1) - x_indent; % unit: nm
if flag==1
    force_R_cor = (handles.data_retract_raw(:,2) - y_indent);
else
    force_R_cor = (handles.data_retract_raw(:,2) - y_indent)*handles.ForceCorFactor; % unit: pN
end
time_R_raw = handles.data_retract_raw(:,3); % unit:s for calculating loading rate dF/dt
handles.data_retract_cor = [extension_R_cor force_R_cor time_R_raw];
handles.data_retract = handles.data_retract_cor; % for smoothing data

axes(handles.axes1); cla reset;
% plot retract curve
handles.h1 = plot(handles.data_retract(:,1),handles.data_retract(:,2),'-r');
hold on
% check if approach file exists: if it does, then plot it
handles.ApproachExist = 0;
if exist(str1, 'file') == 2 %
    handles.ApproachExist = 1;
    handles.data_approach_raw = dlmread(str1);
    extension_A_cor = handles.data_approach_raw(:,1) - x_indent; % unit: nm
    if flag==1
        force_A_cor = (handles.data_approach_raw(:,2) - y_indent); % unit: pN
    else
        force_A_cor = (handles.data_approach_raw(:,2) - y_indent)*handles.ForceCorFactor; % unit: pN
    end
    time_A_raw = handles.data_approach_raw(:,3); % unit:s for calculating loading rate dF/dt
    handles.data_approach_cor = [extension_A_cor force_A_cor time_A_raw];
    handles.data_approach = handles.data_approach_cor;
    
    switch get(get(handles.shift_or_not,'SelectedObject'),'Tag')
        case 'shift_true'
            if ~isempty(handles.shift_dist)
                handles.shiftDist = str2double(get(handles.shift_dist,'String'));
            end
            handles.h2 = plot(handles.data_approach(:,1)+handles.shiftDist, handles.data_approach(:,2),'-b');
        case 'shift_false'
            handles.h2 = plot(handles.data_approach(:,1), handles.data_approach(:,2),'-b');
    end
end
hold off
xlabel('Extension (nm)');ylabel('Force (pN)');
title(['Force vs Extension (# ' num2str(Index_curves(handles.ForceIndex)) ')'])
axis(handles.axisRangeZoomed);

handles.unfolding_num = 0;
handles.refolding_num = 0;
handles.rupture_num = 0;
set(handles.two_curves,'Value',1)
guidata(hObject, handles);



% --- Executes on button press in start_point_1.
function start_point_1_Callback(hObject, eventdata, handles)
% hObject    handle to start_point_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data_retract = handles.data_retract;
dcm_obj = datacursormode(gcf);
info = getCursorInfo(dcm_obj);
p0 = info.Position;
% index0 = data_retract(:,1)==p0(1,1);
index0 = all(bsxfun(@eq,p0,data_retract(:,1:2)),2);
StartPoint = data_retract(index0,:);

handles.StartPoint = StartPoint;
guidata(hObject, handles);


% --- Executes on button press in unfold_point.
function unfold_point_Callback(hObject, eventdata, handles)
% hObject    handle to unfold_point (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles.unfolding_num = handles.unfolding_num + 1;
data_retract = handles.data_retract;
dcm_obj = datacursormode(gcf);
info = getCursorInfo(dcm_obj);
p1 = info.Position;
% index1 = data_retract(:,1)==p1(1,1);
index1 = all(bsxfun(@eq,p1,data_retract(:,1:2)),2);
UnfoldPoint = data_retract(index1,:);
StartPoint = handles.StartPoint;
dF = UnfoldPoint(1,2)-StartPoint(1,2);
dt = UnfoldPoint(1,3)-StartPoint(1,3);
loading_rate = dF/dt;

HozLineX = UnfoldPoint(1,1):0.05:UnfoldPoint(1,1)+100;
HozLineY = ones(size(HozLineX))*UnfoldPoint(1,2);
[xout,yout] = intersections(data_retract(:,1),data_retract(:,2),HozLineX,HozLineY,1);
% average the intersection points, after remobing the 1st intersection
% points, which is the unfolding point
if length(xout)>=3
    xout_mean = mean(xout(2:end-1)); yout_mean = mean(yout(2:end-1));
elseif length(xout)>1
    xout_mean = xout(2); yout_mean = yout(2);
else
    xout_mean = xout(1); yout_mean = yout(1);
end

hold on
plot(HozLineX,HozLineY,':m','LineWidth',1.5);
% plot(xout,yout,'r.','markersize',18)
plot(xout_mean, yout_mean, 'og','markersize',4,'MarkerFaceColor','g')
hold off

set(handles.unfold_force,'string',num2str(UnfoldPoint(1,2),'%5.2f'));
set(handles.loading_rate_1,'string',num2str(loading_rate,'%5.2f'));
set(handles.extension1,'string',num2str(xout_mean-UnfoldPoint(1,1),'%5.2f'));

handles.UnfoldPoint = UnfoldPoint;
handles.IntersectPoint = [xout_mean yout_mean];
handles.UnfoldDistance = xout_mean-UnfoldPoint(1,1);
guidata(hObject, handles);


% --- Executes on button press in intersect_point_1.
function intersect_point_1_Callback(hObject, eventdata, handles)
% hObject    handle to intersect_point_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dcm_obj = datacursormode(gcf);
info = getCursorInfo(dcm_obj);
IntersectPoint = info.Position;
UnfoldPoint = handles.UnfoldPoint;
unfold_distance = IntersectPoint(1,1)-UnfoldPoint(1,1);
set(handles.extension1,'string',num2str(unfold_distance,'%5.2f'));

handles.IntersectPoint = IntersectPoint;
handles.UnfoldDistance = unfold_distance;
guidata(hObject, handles);


% --- Executes on button press in save_refold.
function save_unfold_Callback(hObject, eventdata, handles)
% hObject    handle to save_refold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global filename fid_unfold ExpDate

handles.unfolding_num = handles.unfolding_num+1;
SpeedValue = handles.SpeedValue;
Index_curves = handles.Index_curves;
ForceNum = Index_curves(handles.ForceIndex);
StartPoint = handles.StartPoint;
UnfoldPoint = handles.UnfoldPoint;
unfold_distance = handles.UnfoldDistance;
dF = UnfoldPoint(1,2)-StartPoint(1,2);
dt = UnfoldPoint(1,3)-StartPoint(1,3);
loading_rate = dF/dt;
% dx = UnfoldPoint(1,1)-StartPoint(1,1);

Force_LR_unfold_all = ['data/00 v' num2str(SpeedValue) '_unfold.txt'];
if exist(Force_LR_unfold_all,'file')==2 % if file exists
    fid_unfold = fopen(Force_LR_unfold_all, 'a+'); % append new results
else
    fid_unfold = fopen(Force_LR_unfold_all, 'w+');
    fprintf(fid_unfold, '%9s   %s  %s  %s  %s  %s    %s      %s\r\n', 'LR(N/s)', 'Unfold_F(N)', 'Unfold_D(m)', 'Time(s)', 'Curve#', 'FileName', 'Date', 'F_cor');
end

fprintf(fid_unfold,'%.4e  %.4e  %11.4e  %7.3f  %5d  %9s  %10s   %.4f\r\n',...
    loading_rate*1e-12, UnfoldPoint(1,2)*1e-12, unfold_distance*1e-9, dt, ForceNum, filename, ExpDate, handles.ForceCorFactor);
fclose(fid_unfold);

xLimits = get(gca,'XLim');
yLimits = get(gca,'YLim');

handles.axisRangeZoomed = [xLimits yLimits];
% [extension force time distance loading-rate]
handles.UnfoldPoints(handles.unfolding_num,:) = [UnfoldPoint, unfold_distance, loading_rate];
guidata(hObject, handles);


% --- Executes on button press in start_point_2.
function start_point_2_Callback(hObject, eventdata, handles)
% hObject    handle to start_point_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in refold_point.
function refold_point_Callback(hObject, eventdata, handles)
% hObject    handle to refold_point (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% update retract data if there is a shift


% --- Executes on button press in intersect_point_2.
function intersect_point_2_Callback(hObject, eventdata, handles)
% hObject    handle to intersect_point_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in save_refold.
function save_refold_Callback(hObject, eventdata, handles)
% hObject    handle to save_refold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in start_point_3.
function start_point_3_Callback(hObject, eventdata, handles)
% hObject    handle to start_point_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in rupture_point.
function rupture_point_Callback(hObject, eventdata, handles)
% hObject    handle to rupture_point (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in save_rupture.
function save_rupture_Callback(hObject, eventdata, handles)
% hObject    handle to save_rupture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in save_figure.
function save_figure_Callback(hObject, eventdata, handles)
% hObject    handle to save_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% create a new figure for saving
global filename file_path ExpDate
Index_curves = handles.Index_curves;
ForceNum = Index_curves(handles.ForceIndex);
formatSpec = '%5.2f'; % set display format for parameters
InfoOnFig = {};
fig = figure('visible','off');
set(fig,'units', 'normalized','outerposition',[0.1 0.1 0.8 0.8]);
handles.h1 = plot(handles.data_retract(:,1),handles.data_retract(:,2),'-r');
hold on

if handles.ApproachExist == 1 % if approach file exists
    switch get(get(handles.shift_or_not,'SelectedObject'),'Tag')
        case 'shift_true'
            if ~isempty(handles.shift_dist)
                handles.shiftDist = str2double(get(handles.shift_dist,'String'));
            end
            handles.h2 = plot(handles.data_approach(:,1)+handles.shiftDist, handles.data_approach(:,2),'-b');
        case 'shift_false'
            handles.h2 = plot(handles.data_approach(:,1), handles.data_approach(:,2),'-b');
    end
end

% plot unfolding points
if (handles.unfolding_num>0)
    for i=1:handles.unfolding_num
        hold on
        plot(handles.UnfoldPoints(i,1),handles.UnfoldPoints(i,2),'om','MarkerFaceColor','m','MarkerSize',4);
        InfoOnFig = [InfoOnFig{:},{['F_{u' num2str(i) '}=', num2str(handles.UnfoldPoints(i,2),formatSpec),' pN, ',...
            'D_{u' num2str(i) '}=', num2str(handles.UnfoldPoints(i,4),formatSpec),' nm, ',...
            'LR_{u' num2str(i) '}=', num2str(handles.UnfoldPoints(i,5),formatSpec),' pN/s']}];
    end
end
% plot refolding points
if (handles.refolding_num>0)
    for i=1:handles.refolding_num
        hold on
        plot(handles.RefoldPoints(i,1),handles.RefoldPoints(i,2),'og','MarkerFaceColor','g','MarkerSize',4);
        InfoOnFig = [InfoOnFig{:},{['F_{r' num2str(i) '}=', num2str(handles.RefoldPoints(i,2),formatSpec),' pN, ',...
            'D_{r' num2str(i) '}=', num2str(handles.RefoldPoints(i,4),formatSpec),' nm, ',...
            'LR_{r' num2str(i) '}=', num2str(handles.RefoldPoints(i,5),formatSpec),' pN/s']}];
    end
end
% plot rupture point
if (handles.rupture_num>0)
    hold on
    plot(handles.RupturePoints(:,1),handles.RupturePoint(:,2),'oc','MarkerFaceColor','c','MarkerSize',4);
    InfoOnFig = [InfoOnFig{:},{['F_{ru}=', num2str(handles.RupturePoints(1,2),formatSpec),' pN, ',...
        'D_{ru}=', num2str(handles.RupturePoints(1,4),formatSpec),' nm, ',...
        'LR_{ru}=', num2str(handles.RupturePoints(1,5),formatSpec),' pN/s']}];
end
hold off
xlabel('Extension (nm)');ylabel('Force (pN)');
title(['Force vs Extension (# ' num2str(ForceNum) ')'])
if ~isempty(InfoOnFig)
    annotation('textbox',[0.15 0.7 0.4 0.15],'String',InfoOnFig,'Color','m','BackgroundColor','w')
end
axis(handles.axisRange);
grid off

% warning off MATLAB:MKDIR:DirectoryExists
if ~exist(['data/' ExpDate '/00 Figures'],'dir') % create a new subfolder, if not exist
    mkdir(['data/' ExpDate '/00 Figures']);
end

FigPath = ['data/' ExpDate '/00 Figures/'];
print(gcf,'-dpng',[FigPath filename ' - ' num2str(ForceNum) '.png'])

% save corrected, smoothed data, and delete unsmoothed data file(s)
single_cycle_smooth_file = [file_path 'DFT '];
single_cycle_file = [file_path 'Dist_Force_Time '];
if handles.smoothSpan ~= 1 || handles.ForceCorFactor ~= 1
    str1 = [single_cycle_smooth_file num2str(ForceNum) 'A.txt'];
    str2 = [single_cycle_smooth_file num2str(ForceNum) 'R.txt'];
    if handles.ApproachExist == 1
        fid1 = fopen(str1,'w');
        fid2 = fopen(str2,'w');
        fprintf(fid1,'%9.3f %6.2f %8.3f\r\n',handles.data_approach');
        fprintf(fid2,'%9.3f %6.2f %8.3f\r\n',handles.data_retract');
    else
        fid2 = fopen(str2,'w');
        fprintf(fid2,'%9.3f %6.2f %8.3f\r\n',handles.data_retract');
    end
    
    str3 = [single_cycle_file num2str(ForceNum) 'A.txt'];
    str4 = [single_cycle_file num2str(ForceNum) 'R.txt'];
    if exist(str3, 'file')==2
        delete(str3);
    end
    if exist(str4, 'file')==2
        delete(str4);
    end
end
fclose('all');


function file_name_Callback(hObject, eventdata, handles)
% hObject    handle to file_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filename = get(handles.file_name,'String');
set(handles.speed_value,'String',num2str(filename));

handles.filename = filename;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function file_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to file_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in update_range.
function update_range_Callback(hObject, eventdata, handles)
% hObject    handle to update_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
xLimits(1,1) = str2double(get(handles.x_min,'String'));
xLimits(1,2) = str2double(get(handles.x_max,'String'));
yLimits(1,1) = str2double(get(handles.y_min,'String'));
yLimits(1,2) = str2double(get(handles.y_max,'String'));
handles.axisRange = [xLimits yLimits];
axis(handles.axisRange);
handles.axisRangeZoomed = [xLimits yLimits];

guidata(hObject, handles);


function x_min_Callback(hObject, eventdata, handles)
% hObject    handle to x_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
xLimits = get(gca,'XLim'); yLimits = get(gca,'YLim');
xLimits(1,1) = str2double(get(handles.x_min,'String'));
xLimits(1,2) = xLimits(1,1)+handles.XSpan;
set(handles.x_max,'string',num2str(xLimits(2),'%5d'));
handles.axisRange = [xLimits yLimits];
axis(handles.axisRange);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function x_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function x_max_Callback(hObject, eventdata, handles)
% hObject    handle to x_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
xLimits = get(gca,'XLim'); yLimits = get(gca,'YLim');
xLimits(1,2) = str2double(get(handles.x_max,'String'));
xLimits(1,1) = xLimits(1,2)-handles.XSpan;
set(handles.x_min,'string',num2str(xLimits(1),'%5d'));
handles.axisRange = [xLimits yLimits];
axis(handles.axisRange);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function x_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function y_min_Callback(hObject, eventdata, handles)
% hObject    handle to y_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
xLimits = get(gca,'XLim'); yLimits = get(gca,'YLim');
yLimits(1,1) = str2double(get(handles.y_min,'String'));
yLimits(1,2) = yLimits(1,1)+handles.YSpan;
set(handles.y_max,'string',num2str(yLimits(2),'%5d'));
handles.axisRange = [xLimits yLimits];
axis(handles.axisRange);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function y_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function y_max_Callback(hObject, eventdata, handles)
% hObject    handle to y_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
xLimits = get(gca,'XLim'); yLimits = get(gca,'YLim');
yLimits(1,2) = str2double(get(handles.y_max,'String'));
yLimits(1,1) = yLimits(1,2)-handles.YSpan;
set(handles.y_min,'string',num2str(yLimits(1),'%5d'));
handles.axisRange = [xLimits yLimits];
axis(handles.axisRange);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function y_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function x_span_Callback(hObject, eventdata, handles)
% hObject    handle to x_span (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
xLimits = get(gca,'XLim'); yLimits = get(gca,'YLim');
handles.XSpan = str2double(get(handles.x_span,'String'));
xLimits(1,2) = xLimits(1,1)+handles.XSpan;
set(handles.x_max,'string',num2str(xLimits(2),'%5d'));
handles.axisRange = [xLimits yLimits];
axis(handles.axisRange);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function x_span_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_span (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function y_span_Callback(hObject, eventdata, handles)
% hObject    handle to y_span (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
xLimits = get(gca,'XLim'); yLimits = get(gca,'YLim');
handles.YSpan = str2double(get(handles.y_span,'String'));
yLimits(1,2) = yLimits(1,1)+handles.YSpan;
set(handles.y_max,'string',num2str(yLimits(2),'%5d'));
handles.axisRange = [xLimits yLimits];
axis(handles.axisRange);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function y_span_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_span (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function unfold_force_Callback(hObject, eventdata, handles)
% hObject    handle to unfold_force (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function unfold_force_CreateFcn(hObject, eventdata, handles)
% hObject    handle to unfold_force (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function refold_force_Callback(hObject, eventdata, handles)
% hObject    handle to refold_force (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function refold_force_CreateFcn(hObject, eventdata, handles)
% hObject    handle to refold_force (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function rupture_force_Callback(hObject, eventdata, handles)
% hObject    handle to rupture_force (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function rupture_force_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rupture_force (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function loading_rate_1_Callback(hObject, eventdata, handles)
% hObject    handle to loading_rate_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function loading_rate_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loading_rate_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function loading_rate_2_Callback(hObject, eventdata, handles)
% hObject    handle to loading_rate_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function loading_rate_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loading_rate_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function loading_rate_3_Callback(hObject, eventdata, handles)
% hObject    handle to loading_rate_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function loading_rate_3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loading_rate_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function extension1_Callback(hObject, eventdata, handles)
% hObject    handle to extension1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function extension1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to extension1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function extension2_Callback(hObject, eventdata, handles)
% hObject    handle to extension2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function extension2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to extension2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function extension3_Callback(hObject, eventdata, handles)
% hObject    handle to extension3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function extension3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to extension3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in previous_file.
function previous_file_Callback(hObject, eventdata, handles)
% hObject    handle to previous_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global file_path Index_curves
if handles.ForceIndex>1 && length(Index_curves)>1
    % read data file(s)
    handles.ForceIndex = handles.ForceIndex - 1;
else
    handles.ForceIndex = handles.ForceIndex;
end
ForceNum = Index_curves(handles.ForceIndex);

single_cycle_file = [file_path 'Dist_Force_Time '];
single_cycle_smooth_file = [file_path 'DFT '];

str1 = [single_cycle_file num2str(ForceNum) 'A.txt'];
str2 = [single_cycle_file num2str(ForceNum) 'R.txt'];

% check if the data is smoothed or not
flag = 0;
if exist(str2, 'file') ~= 2
    str1 = [single_cycle_smooth_file num2str(ForceNum) 'A.txt'];
    str2 = [single_cycle_smooth_file num2str(ForceNum) 'R.txt'];
    flag = 1;
end
handles.data_retract_raw = [];
handles.data_approach_raw = [];
handles.data_retract_cor = [];
handles.data_approach_cor = [];
handles.data_retract = [];
handles.data_approach = [];

x_indent = 0;
y_indent = 0;

handles.data_retract_raw = dlmread(str2);
extension_R_cor = handles.data_retract_raw(:,1) - x_indent; % unit: nm
if flag==1
    force_R_cor = (handles.data_retract_raw(:,2) - y_indent);
else
    force_R_cor = (handles.data_retract_raw(:,2) - y_indent)*handles.ForceCorFactor; % unit: pN
end
time_R_raw = handles.data_retract_raw(:,3); % unit:s for calculating loading rate dF/dt
handles.data_retract_cor = [extension_R_cor force_R_cor time_R_raw];
handles.data_retract = handles.data_retract_cor; % for smoothing data

axes(handles.axes1); cla reset;
% plot retract curve
handles.h1 = plot(handles.data_retract(:,1),handles.data_retract(:,2),'-r');
hold on
% check if approach file exists: if it does, then plot it
handles.ApproachExist = 0;
if exist(str1, 'file') == 2 %
    handles.ApproachExist = 1;
    handles.data_approach_raw = dlmread(str1);
    extension_A_cor = handles.data_approach_raw(:,1) - x_indent; % unit: nm
    if flag==1
        force_A_cor = (handles.data_approach_raw(:,2) - y_indent); % unit: pN
    else
        force_A_cor = (handles.data_approach_raw(:,2) - y_indent)*handles.ForceCorFactor; % unit: pN
    end
    time_A_raw = handles.data_approach_raw(:,3); % unit:s for calculating loading rate dF/dt
    handles.data_approach_cor = [extension_A_cor force_A_cor time_A_raw];
    handles.data_approach = handles.data_approach_cor;
    
    switch get(get(handles.shift_or_not,'SelectedObject'),'Tag')
        case 'shift_true'
            if ~isempty(handles.shift_dist)
                handles.shiftDist = str2double(get(handles.shift_dist,'String'));
            end
            handles.h2 = plot(handles.data_approach(:,1)+handles.shiftDist, handles.data_approach(:,2),'-b');
        case 'shift_false'
            handles.h2 = plot(handles.data_approach(:,1), handles.data_approach(:,2),'-b');
    end

end
hold off
xlabel('Extension (nm)');ylabel('Force (pN)');
title(['Force vs Extension (# ' num2str(Index_curves(handles.ForceIndex)) ')'])
axis(handles.axisRangeZoomed);
%     grid on
%     datacursormode on;

handles.unfolding_num = 0;
handles.refolding_num = 0;
handles.rupture_num = 0;
guidata(hObject, handles);


% --- Executes when selected object is changed in curve_selection.
function curve_selection_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in curve_selection
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
Index_curves = handles.Index_curves;
xLimits = get(gca,'XLim'); yLimits = get(gca,'YLim');

handles.h1 = plot(handles.data_retract(:,1),handles.data_retract(:,2),'-r');
hold on
switch get(get(handles.shift_or_not,'SelectedObject'),'Tag')
    case 'shift_true'
        if ~isempty(handles.shift_dist)
            handles.shiftDist = str2double(get(handles.shift_dist,'String'));
        end
        handles.h2 = plot(handles.data_approach(:,1)+handles.shiftDist, handles.data_approach(:,2),'-b');
    case 'shift_false'
        handles.h2 = plot(handles.data_approach(:,1), handles.data_approach(:,2),'-b');
end
hold off

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'two_curves'

    case 'retract_only'
        delete(handles.h2)
    case 'approach_only'
        delete(handles.h1)
end
xlabel('Extension (nm)');ylabel('Force (pN)');
title(['Force vs Extension (# ' num2str(Index_curves(handles.ForceIndex)) ')'])
axis([xLimits yLimits]);


function go_curve_Callback(hObject, eventdata, handles)
% hObject    handle to go_curve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global file_path
single_cycle_file = [ file_path 'Dist_Force_Time ']; % one approaching/retracting process
single_cycle_smooth_file = [ file_path 'DFT ']; % one approaching/retracting process

Index_curves = handles.Index_curves;
ForceNum = str2double(get(handles.go_curve,'String'));
ForceIndex = find(Index_curves==ForceNum);
handles.ForceIndex = ForceIndex;

str1 = [single_cycle_file num2str(ForceNum) 'A.txt'];
str2 = [single_cycle_file num2str(ForceNum) 'R.txt'];

% check if the data is smoothed or not
flag = 0;
if exist(str2, 'file') ~= 2
    str1 = [single_cycle_smooth_file num2str(ForceNum) 'A.txt'];
    str2 = [single_cycle_smooth_file num2str(ForceNum) 'R.txt'];
    flag = 1;
end

handles.data_retract_raw = [];
handles.data_approach_raw = [];
handles.data_retract_cor = [];
handles.data_approach_cor = [];
handles.data_retract = [];
handles.data_approach = [];

x_indent = 0;
y_indent = 0;

handles.data_retract_raw = dlmread(str2);
extension_R_cor = handles.data_retract_raw(:,1) - x_indent; % unit: nm
if flag==1
    force_R_cor = (handles.data_retract_raw(:,2) - y_indent);
else
    force_R_cor = (handles.data_retract_raw(:,2) - y_indent)*handles.ForceCorFactor; % unit: pN
end
time_R_raw = handles.data_retract_raw(:,3); % unit:s for calculating loading rate dF/dt
handles.data_retract_cor = [extension_R_cor force_R_cor time_R_raw];
handles.data_retract = handles.data_retract_cor; % for smoothing data

axes(handles.axes1); cla reset;
% plot retract curve
plot(handles.data_retract_cor(:,1),handles.data_retract_cor(:,2),'-r')

% check if approach file exists: if it does, then plot it
handles.ApproachExist = 0;
if exist(str1, 'file') == 2 %
    handles.ApproachExist = 1;
    handles.data_approach_raw = dlmread(str1);
    extension_A_cor = handles.data_approach_raw(:,1) - x_indent; % unit: nm
    if flag==1
        force_A_cor = (handles.data_approach_raw(:,2) - y_indent); % unit: pN
    else
        force_A_cor = (handles.data_approach_raw(:,2) - y_indent)*handles.ForceCorFactor; % unit: pN
    end
    time_A_raw = handles.data_approach_raw(:,3); % unit:s for calculating loading rate dF/dt
    handles.data_approach_cor = [extension_A_cor force_A_cor time_A_raw];
    handles.data_approach = handles.data_approach_cor;
    
    hold on
    if handles.shift_true==1
        handles.h2 = plot(handles.data_approach_cor(:,1)+handles.shiftDist, handles.data_approach_cor(:,2),'-b');
    else
        handles.h2 = plot(handles.data_approach_cor(:,1), handles.data_approach_cor(:,2),'-b');
    end
    hold off
end
xlabel('Extension (nm)');ylabel('Force (pN)');
title(['Force vs Extension (# ' num2str(Index_curves(handles.ForceIndex)) ')'])
axis(handles.axisRange);
% grid on

handles.unfolding_num = 0;
handles.refolding_num = 0;
handles.rupture_num = 0;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function go_curve_CreateFcn(hObject, eventdata, handles)
% hObject    handle to go_curve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function smooth_span_Callback(hObject, eventdata, handles)
% hObject    handle to smooth_span (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Index_curves = handles.Index_curves;
smoothSpan = str2double(get(hObject,'String'));
retract_dft = handles.data_retract_cor;
retract_f_smooth = smooth(retract_dft(:,2), smoothSpan);

xLimits = get(gca,'XLim'); yLimits = get(gca,'YLim');
axes(handles.axes1); cla reset;
if handles.ApproachExist == 1 % if approach file exists
    approach_dft = handles.data_approach_cor;
    switch get(get(handles.shift_or_not,'SelectedObject'),'Tag')
        case 'shift_true'
            if ~isempty(handles.shift_dist)
                handles.shiftDist = str2double(get(handles.shift_dist,'String'));
            end
            approach_dft(:,1) = approach_dft(:,1)+handles.shiftDist;
        case 'shift_false'
            
    end
    approach_f_smooth = smooth(approach_dft(:,2), smoothSpan);
    plot(retract_dft(:,1),retract_f_smooth,'-r',approach_dft(:,1), approach_f_smooth,'-b')
    %     legend('Retract','Approach','Location','NorthEast')
    handles.data_approach(:,2) = approach_f_smooth;
else
    plot(retract_dft(:,1),retract_f_smooth,'-r')
    %     legend('Retract','Location','NorthEast')
end
xlabel('Extension (nm)');ylabel('Force (pN)');
title(['Force vs Extension (# ' num2str(Index_curves(handles.ForceIndex)) ')'])
axis([xLimits yLimits]);

handles.smoothSpan = smoothSpan;
handles.data_retract(:,2) = retract_f_smooth;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function smooth_span_CreateFcn(hObject, eventdata, handles)
% hObject    handle to smooth_span (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in update_force.
function update_force_Callback(hObject, eventdata, handles)
% hObject    handle to update_force (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global filename Index_curves ExpDate
% filename = evalin('base', 'filename');
% Index_curves = handles.Index_curves;
newForceCorFactor = str2double(get(handles.force_cor_factor,'String'));

retract_dft_smooth = handles.data_retract;
retract_f_smooth = retract_dft_smooth(:,2)*newForceCorFactor;
handles.data_retract(:,2) = retract_f_smooth;
% update corrected force for further smoothing
handles.data_retract_cor(:,2) = handles.data_retract_cor(:,2)*newForceCorFactor;

if handles.ApproachExist == 1
    approach_dft_smooth = handles.data_approach;
    approach_f_smooth = approach_dft_smooth(:,2)*newForceCorFactor;
    handles.data_approach(:,2) = approach_f_smooth;
    handles.data_approach_cor(:,2) = handles.data_approach_cor(:,2)*newForceCorFactor;
end

axes(handles.axes1);
% plot retract curve
handles.h1 = plot(retract_dft_smooth(:,1),retract_f_smooth,'-r');
hold on
% check if approach file exists: if it does, then plot it
if handles.ApproachExist == 1
    switch get(get(handles.shift_or_not,'SelectedObject'),'Tag')
        case 'shift_true'
            if ~isempty(handles.shift_dist)
                handles.shiftDist = str2double(get(handles.shift_dist,'String'));
            end
            handles.h2 = plot(approach_dft_smooth(:,1)+handles.shiftDist, approach_f_smooth,'-b');
        case 'shift_false'
            handles.h2 = plot(approach_dft_smooth(:,1), approach_f_smooth,'-b');
    end
end
hold off
xlabel('Extension (nm)');ylabel('Force (pN)');
title(['Force vs Extension (# ' num2str(Index_curves(handles.ForceIndex)) ')'])
axis(handles.axisRange);
% grid on

% save force correct factor for next file
handles.ForceCorFactor = handles.ForceCorFactor*newForceCorFactor;
force_factor_file = ['data/' ExpDate '/Force_Cor_Factor.txt']; % [dataset factor]
fid0 = fopen(force_factor_file,'a+');
fprintf(fid0,'%s %.4f\r\n',filename,handles.ForceCorFactor);
fclose(fid0);
guidata(hObject, handles);

function current_force_Callback(hObject, eventdata, handles)
% hObject    handle to current_force (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function current_force_CreateFcn(hObject, eventdata, handles)
% hObject    handle to current_force (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function force_cor_factor_Callback(hObject, eventdata, handles)
% hObject    handle to force_cor_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function force_cor_factor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to force_cor_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function target_force_Callback(hObject, eventdata, handles)
% hObject    handle to target_force (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of target_force as text
%        str2double(get(hObject,'String')) returns contents of target_force as a double


% --- Executes during object creation, after setting all properties.
function target_force_CreateFcn(hObject, eventdata, handles)
% hObject    handle to target_force (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in get_point.
function get_point_Callback(hObject, eventdata, handles)
% hObject    handle to get_point (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dcm_obj = datacursormode(gcf);
info = getCursorInfo(dcm_obj);
CurrentForce = info.Position(2);

set(handles.current_force,'string',num2str(CurrentForce,'%.2f'));
set(handles.force_cor_factor,'string',num2str(handles.TargetForce/CurrentForce,'%.4f'));

guidata(hObject, handles);


function del_bad_curve_Callback(hObject, eventdata, handles)
global file_path Index_curves
Index_curves = handles.Index_curves;
ForceNum = Index_curves(handles.ForceIndex);
str1 = [file_path 'Dist_Force_Time ' num2str(ForceNum) 'A.txt'];
str2 = [file_path 'Dist_Force_Time ' num2str(ForceNum) 'R.txt'];
if exist(str1, 'file')==2
    delete(str1);
end
if exist(str2, 'file')==2
    delete(str2);
end
% update Index_curves, by removing the ForceNum just deleted
% Index_curves(handles.ForceIndex) = [];
handles.Index_curves(handles.Index_curves == ForceNum) = [];
handles.ForceIndex = handles.ForceIndex-1;
assignin('base', 'Index_curves', handles.Index_curves)
guidata(hObject, handles);


function set_average_points_Callback(hObject, eventdata, handles)
Index_curves = handles.Index_curves;
smoothSpan = handles.smoothSpan;
retract_dft = handles.data_retract_cor;
retract_f_smooth = smooth(retract_dft(:,2), smoothSpan);

xLimits = get(gca,'XLim'); yLimits = get(gca,'YLim');
axes(handles.axes1); cla reset;
handles.h1 = plot(retract_dft(:,1),retract_f_smooth,'-r');
hold on
if handles.ApproachExist == 1 % if approach file exists
    approach_dft = handles.data_approach_cor;
    approach_f_smooth = smooth(approach_dft(:,2), smoothSpan);
    switch get(get(handles.shift_or_not,'SelectedObject'),'Tag')
        case 'shift_true'
            if ~isempty(handles.shift_dist)
                handles.shiftDist = str2double(get(handles.shift_dist,'String'));
            end
            handles.h2 = plot(approach_dft(:,1)+handles.shiftDist, approach_f_smooth,'-b');
        case 'shift_false'
            handles.h2 = plot(approach_dft(:,1), approach_f_smooth,'-b');
    end
    handles.data_approach(:,2) = approach_f_smooth;
else
    plot(retract_dft(:,1),retract_f_smooth,'-r')
end
hold off
xlabel('Extension (nm)');ylabel('Force (pN)');
title(['Force vs Extension (# ' num2str(Index_curves(handles.ForceIndex)) ')'])
axis([xLimits yLimits]);

handles.smoothSpan = smoothSpan;
handles.data_retract(:,2) = retract_f_smooth;
guidata(hObject, handles);


function shift_or_not_SelectionChangedFcn(hObject, eventdata, handles)
Index_curves = handles.Index_curves;
xLimits = get(gca,'XLim'); yLimits = get(gca,'YLim');
handles.shiftDist = str2double(get(handles.shift_dist,'String'));
handles.h1 = plot(handles.data_retract(:,1),handles.data_retract(:,2),'-r');
hold on
switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'shift_true'
        handles.h2 = plot(handles.data_approach(:,1)+handles.shiftDist, handles.data_approach(:,2),'-b');
    case 'shift_false'
        delete(handles.h2)
        handles.h2 = plot(handles.data_approach(:,1), handles.data_approach(:,2),'-b');
end
hold off
xlabel('Extension (nm)');ylabel('Force (pN)');
title(['Force vs Extension (# ' num2str(Index_curves(handles.ForceIndex)) ')'])
axis([xLimits yLimits]);

function shift_or_not_CreateFcn(hObject, eventdata, handles)
% hObject    handle to target_force (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function shift_dist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to unfold_force (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function close_gui_Callback(hObject, eventdata, handles)
close all;


function del_next_Callback(hObject, eventdata, handles)
global file_path
Index_curves = handles.Index_curves;
ForceNum = Index_curves(handles.ForceIndex);
str1 = [file_path 'Dist_Force_Time ' num2str(ForceNum) 'A.txt'];
str2 = [file_path 'Dist_Force_Time ' num2str(ForceNum) 'R.txt'];
if exist(str1, 'file')==2
    delete(str1);
end
if exist(str2, 'file')==2
    delete(str2);
end
% update Index_curves, by removing the ForceNum just deleted
% Index_curves(handles.ForceIndex) = [];
handles.Index_curves(handles.Index_curves == ForceNum) = [];
% Index_curves(handles.Index_curves == ForceNum) = [];
assignin('base', 'Index_curves', handles.Index_curves)

Index_curves = handles.Index_curves;
% if handles.ForceIndex<length(Index_curves) && length(Index_curves)>1
%     handles.ForceIndex = handles.ForceIndex + 1;
% else
%     handles.ForceIndex = handles.ForceIndex;
% end
ForceNum = Index_curves(handles.ForceIndex);

single_cycle_file = [file_path 'Dist_Force_Time '];
single_cycle_smooth_file = [file_path 'DFT '];

str1 = [single_cycle_file num2str(ForceNum) 'A.txt'];
str2 = [single_cycle_file num2str(ForceNum) 'R.txt'];

% check if the data is smoothed or not
flag = 0; % check whether to correct force
if exist(str2, 'file') ~= 2
    str1 = [single_cycle_smooth_file num2str(ForceNum) 'A.txt'];
    str2 = [single_cycle_smooth_file num2str(ForceNum) 'R.txt'];
    flag = 1;
end

handles.data_retract_raw = [];
handles.data_approach_raw = [];
handles.data_retract_cor = [];
handles.data_approach_cor = [];
handles.data_retract = [];
handles.data_approach = [];

x_indent = 0;
y_indent = 0;

handles.data_retract_raw = dlmread(str2);
extension_R_cor = handles.data_retract_raw(:,1) - x_indent; % unit: nm
if flag==1
    force_R_cor = (handles.data_retract_raw(:,2) - y_indent);
else
    force_R_cor = (handles.data_retract_raw(:,2) - y_indent)*handles.ForceCorFactor; % unit: pN
end
time_R_raw = handles.data_retract_raw(:,3); % unit:s for calculating loading rate dF/dt
handles.data_retract_cor = [extension_R_cor force_R_cor time_R_raw];
handles.data_retract = handles.data_retract_cor; % for smoothing data

axes(handles.axes1); cla reset;
% plot retract curve
handles.h1 = plot(handles.data_retract(:,1),handles.data_retract(:,2),'-r');
hold on
% check if approach file exists: if it does, then plot it
handles.ApproachExist = 0;
if exist(str1, 'file') == 2 %
    handles.ApproachExist = 1;
    handles.data_approach_raw = dlmread(str1);
    extension_A_cor = handles.data_approach_raw(:,1) - x_indent; % unit: nm
    if flag==1
        force_A_cor = (handles.data_approach_raw(:,2) - y_indent); % unit: pN
    else
        force_A_cor = (handles.data_approach_raw(:,2) - y_indent)*handles.ForceCorFactor; % unit: pN
    end
    time_A_raw = handles.data_approach_raw(:,3); % unit:s for calculating loading rate dF/dt
    handles.data_approach_cor = [extension_A_cor force_A_cor time_A_raw];
    handles.data_approach = handles.data_approach_cor;
    
    switch get(get(handles.shift_or_not,'SelectedObject'),'Tag')
        case 'shift_true'
            if ~isempty(handles.shift_dist)
                handles.shiftDist = str2double(get(handles.shift_dist,'String'));
            end
            handles.h2 = plot(handles.data_approach(:,1)+handles.shiftDist, handles.data_approach(:,2),'-b');
        case 'shift_false'
            handles.h2 = plot(handles.data_approach(:,1), handles.data_approach(:,2),'-b');
    end
end
hold off
xlabel('Extension (nm)');ylabel('Force (pN)');
title(['Force vs Extension (# ' num2str(Index_curves(handles.ForceIndex)) ')'])
axis(handles.axisRangeZoomed);

handles.unfolding_num = 0;
handles.refolding_num = 0;
handles.rupture_num = 0;
set(handles.two_curves,'Value',1)


% del_bad_curve_Callback(handles.del_bad_curve, eventdata, handles)

% next_file_Callback(handles.next_file, eventdata, handles)
% handles.ForceIndex = handles.ForceIndex-1;
guidata(hObject, handles);


function go_first_Callback(hObject, eventdata, handles)
global file_path
single_cycle_file = [ file_path 'Dist_Force_Time ']; % one approaching/retracting process
single_cycle_smooth_file = [ file_path 'DFT ']; % one approaching/retracting process

Index_curves = handles.Index_curves;
ForceNum = Index_curves(1);
ForceIndex = 1;
handles.ForceIndex = ForceIndex;

str1 = [single_cycle_file num2str(ForceNum) 'A.txt'];
str2 = [single_cycle_file num2str(ForceNum) 'R.txt'];

% check if the data is smoothed or not
flag = 0;
if exist(str2, 'file') ~= 2
    str1 = [single_cycle_smooth_file num2str(ForceNum) 'A.txt'];
    str2 = [single_cycle_smooth_file num2str(ForceNum) 'R.txt'];
    flag = 1;
end

handles.data_retract_raw = [];
handles.data_approach_raw = [];
handles.data_retract_cor = [];
handles.data_approach_cor = [];
handles.data_retract = [];
handles.data_approach = [];

x_indent = 0;
y_indent = 0;

handles.data_retract_raw = dlmread(str2);
extension_R_cor = handles.data_retract_raw(:,1) - x_indent; % unit: nm
if flag==1
    force_R_cor = (handles.data_retract_raw(:,2) - y_indent);
else
    force_R_cor = (handles.data_retract_raw(:,2) - y_indent)*handles.ForceCorFactor; % unit: pN
end
time_R_raw = handles.data_retract_raw(:,3); % unit:s for calculating loading rate dF/dt
handles.data_retract_cor = [extension_R_cor force_R_cor time_R_raw];
handles.data_retract = handles.data_retract_cor; % for smoothing data

axes(handles.axes1); cla reset;
% plot retract curve
plot(handles.data_retract_cor(:,1),handles.data_retract_cor(:,2),'-r')

% check if approach file exists: if it does, then plot it
handles.ApproachExist = 0;
if exist(str1, 'file') == 2 %
    handles.ApproachExist = 1;
    handles.data_approach_raw = dlmread(str1);
    extension_A_cor = handles.data_approach_raw(:,1) - x_indent; % unit: nm
    if flag==1
        force_A_cor = (handles.data_approach_raw(:,2) - y_indent); % unit: pN
    else
        force_A_cor = (handles.data_approach_raw(:,2) - y_indent)*handles.ForceCorFactor; % unit: pN
    end
    time_A_raw = handles.data_approach_raw(:,3); % unit:s for calculating loading rate dF/dt
    handles.data_approach_cor = [extension_A_cor force_A_cor time_A_raw];
    handles.data_approach = handles.data_approach_cor;
    
    hold on
    if handles.shift_true==1
        handles.h2 = plot(handles.data_approach_cor(:,1)+handles.shiftDist, handles.data_approach_cor(:,2),'-b');
    else
        handles.h2 = plot(handles.data_approach_cor(:,1), handles.data_approach_cor(:,2),'-b');
    end
    hold off
end
xlabel('Extension (nm)');ylabel('Force (pN)');
title(['Force vs Extension (# ' num2str(Index_curves(handles.ForceIndex)) ')'])
axis(handles.axisRange);
% grid on

handles.unfolding_num = 0;
handles.refolding_num = 0;
handles.rupture_num = 0;
guidata(hObject, handles);

function go_last_Callback(hObject, eventdata, handles)
global file_path
single_cycle_file = [ file_path 'Dist_Force_Time ']; % one approaching/retracting process
single_cycle_smooth_file = [ file_path 'DFT ']; % one approaching/retracting process

Index_curves = handles.Index_curves;
ForceNum = Index_curves(end);
ForceIndex = length(Index_curves);
handles.ForceIndex = ForceIndex;

str1 = [single_cycle_file num2str(ForceNum) 'A.txt'];
str2 = [single_cycle_file num2str(ForceNum) 'R.txt'];

% check if the data is smoothed or not
flag = 0;
if exist(str2, 'file') ~= 2
    str1 = [single_cycle_smooth_file num2str(ForceNum) 'A.txt'];
    str2 = [single_cycle_smooth_file num2str(ForceNum) 'R.txt'];
    flag = 1;
end

handles.data_retract_raw = [];
handles.data_approach_raw = [];
handles.data_retract_cor = [];
handles.data_approach_cor = [];
handles.data_retract = [];
handles.data_approach = [];

x_indent = 0;
y_indent = 0;

handles.data_retract_raw = dlmread(str2);
extension_R_cor = handles.data_retract_raw(:,1) - x_indent; % unit: nm
if flag==1
    force_R_cor = (handles.data_retract_raw(:,2) - y_indent);
else
    force_R_cor = (handles.data_retract_raw(:,2) - y_indent)*handles.ForceCorFactor; % unit: pN
end
time_R_raw = handles.data_retract_raw(:,3); % unit:s for calculating loading rate dF/dt
handles.data_retract_cor = [extension_R_cor force_R_cor time_R_raw];
handles.data_retract = handles.data_retract_cor; % for smoothing data

axes(handles.axes1); cla reset;
% plot retract curve
plot(handles.data_retract_cor(:,1),handles.data_retract_cor(:,2),'-r')

% check if approach file exists: if it does, then plot it
handles.ApproachExist = 0;
if exist(str1, 'file') == 2 %
    handles.ApproachExist = 1;
    handles.data_approach_raw = dlmread(str1);
    extension_A_cor = handles.data_approach_raw(:,1) - x_indent; % unit: nm
    if flag==1
        force_A_cor = (handles.data_approach_raw(:,2) - y_indent); % unit: pN
    else
        force_A_cor = (handles.data_approach_raw(:,2) - y_indent)*handles.ForceCorFactor; % unit: pN
    end
    time_A_raw = handles.data_approach_raw(:,3); % unit:s for calculating loading rate dF/dt
    handles.data_approach_cor = [extension_A_cor force_A_cor time_A_raw];
    handles.data_approach = handles.data_approach_cor;
    
    hold on
    if handles.shift_true==1
        handles.h2 = plot(handles.data_approach_cor(:,1)+handles.shiftDist, handles.data_approach_cor(:,2),'-b');
    else
        handles.h2 = plot(handles.data_approach_cor(:,1), handles.data_approach_cor(:,2),'-b');
    end
    hold off
end
xlabel('Extension (nm)');ylabel('Force (pN)');
title(['Force vs Extension (# ' num2str(Index_curves(handles.ForceIndex)) ')'])
axis(handles.axisRange);
% grid on

handles.unfolding_num = 0;
handles.refolding_num = 0;
handles.rupture_num = 0;
guidata(hObject, handles);
