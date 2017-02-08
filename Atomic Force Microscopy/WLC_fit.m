% Nonlinear fitting about force-extension curve by WLC model
% Wei Zhang (wez311@lehigh.edu)
% Lehigh University
% units: Lc: nm; Lp: nm; x: nm; F: pN; spring: N/m; slope: m/V

function varargout = WLC_fit(varargin)
% WLC_FIT MATLAB code for WLC_fit.fig
%      WLC_FIT, by itself, creates a new WLC_FIT or raises the existing
%      singleton*.
%
%      H = WLC_FIT returns the handle to a new WLC_FIT or the handle to
%      the existing singleton*.
%
%      WLC_FIT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WLC_FIT.M with the given input arguments.
%
%      WLC_FIT('Property','Value',...) creates a new WLC_FIT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before WLC_fit_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      update.  All inputs are passed to WLC_fit_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help WLC_fit

% Last Modified by GUIDE v2.5 18-Jan-2017 14:37:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @WLC_fit_OpeningFcn, ...
    'gui_OutputFcn',  @WLC_fit_OutputFcn, ...
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


% --- Executes just before WLC_fit is made visible.
function WLC_fit_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to WLC_fit (see VARARGIN)

clc;
global kB T;
kB = 1.3806503*1e-23; % Boltzmann constant
T = 273.15+24; % absolute temperature

% Choose default command line output for WLC_fit
handles.output = hObject;

% initilize parameters
handles.PeakNum = 1;
handles.unit_K = 'N/m';
handles.unit_s = 'm/V';
set(handles.unit_spring,'string',handles.unit_K);
set(handles.unit_slope,'string',handles.unit_s);

handles.spring = 0.01;
handles.slope = 5e-8;
set(handles.spring_value,'string',num2str(handles.spring));
set(handles.slope_value,'string',num2str(handles.slope));

% Initilize linear fit parameters (baseline) and origin
handles.p = [0 0]; % linear fit parameters
handles.origin = [0 0];
%% load list box
dir_path = pwd;
source_files = [dir(fullfile(dir_path, '*.xlsx')); dir(fullfile(dir_path, '*.xls'))];
for i =1:length(source_files)
    handles.excellist{i} = source_files(i).name; %lists each of them in a cell
end
set(handles.FileList,'String',handles.excellist,'Value',1)

handles.fit_max_force = 400; % unit: pN
set(handles.y_max_fit,'string',num2str(handles.fit_max_force));

handles.source_files = source_files;
handles.file_path = 'data/';
handles.fig_path = 'figure/';
if ~exist(handles.file_path, 'dir')
    mkdir(fn)
end
if ~exist(handles.fig_path, 'dir')
    mkdir(handles.fig_path)
end
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = WLC_fit_OutputFcn(hObject, eventdata, handles)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in FileList.
function FileList_Callback(hObject, eventdata, handles)
% hObject    handle to FileList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns FileList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FileList
global FilePath CsvFile FileName Index_curves ForceIndex slope spring;
FilePath = handles.file_path;
cla reset % reset figure
handles.PeakNum = 1;
if isfield(handles, 'Lc_Lp')
    handles = rmfield(handles,'Lc_Lp');
end

list = get(handles.FileList,'string'); %get the picture list
handles.FileNum = get(handles.FileList,'value'); % get which one is selected.
handles.File = fullfile(pwd,cell2mat(list(handles.FileNum)));
[~, FileName, ~] = fileparts(handles.File);
set(handles.filename,'string',FileName);

% read raw data
[handles.data_raw, header, ~] = xlsread(FileName, 1); % header is [force, indent, force, indent, ...]
Index_curves = zeros(1, length(header)/2);
for i=1:length(header)/2
    Index_curves(i) = str2double(header{2*i}(7:9));
end

ForceIndex = 1;
spring = handles.spring; % N/m
slope = handles.slope; % m/V

% AFM data
CsvFile = [FilePath 'Force' num2str(Index_curves(ForceIndex),'%03d') '.csv']; % [Extension Force] (nm, pN)
if exist(CsvFile,'file')
    data = dlmread(CsvFile,',',1,0);
    extension = data(:,1);
    force = data(:,2);
else
    data_retract_raw = handles.data_raw(:,2*(ForceIndex-1)+2:-1:2*(ForceIndex-1)+1); % [Extension Force]
    extension_raw = data_retract_raw(:,1); % unit: m
    force = data_retract_raw(:,2)*spring*slope*1e12; % unit: pN
    extension = extension_raw*1e9 - force/spring*1e-3; % unit: nm, consider the cantilever shape change
end
% plot the retract curve
datacursormode on
plot(extension, force)
xlabel('Extension (nm)');ylabel('Force (pN)');
title(['Force vs Extension (# ' num2str(Index_curves(ForceIndex)) ')'])
ylim([-150 max(force)+50])

xLimits = get(gca,'XLim'); yLimits = get(gca,'YLim');
XSpan = xLimits(2)-xLimits(1);
YSpan = yLimits(2)-yLimits(1);

set(handles.x_min,'string',num2str(xLimits(1),'%.0f'));
set(handles.x_max,'string',num2str(xLimits(2),'%.0f'));
set(handles.y_min,'string',num2str(yLimits(1),'%.0f'));
set(handles.y_max,'string',num2str(yLimits(2),'%.0f'));
set(handles.x_span,'string',num2str(XSpan,'%.0f'));
set(handles.y_span,'string',num2str(YSpan,'%.0f'));

handles.XSpan = XSpan; handles.YSpan = YSpan;
handles.axisRange = [xLimits yLimits];
handles.axisRangeZoomed = [xLimits yLimits];
handles.data_retract_cor = [extension force];
handles.data_retract = handles.data_retract_cor; % for smoothing data
handles.Index_curves = Index_curves;
handles.ForceIndex = ForceIndex;
handles.origin = zeros(1,2,length(Index_curves)); % save origin for each curve
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function FileList_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function next_curve_Callback(hObject, eventdata, handles)
global FilePath CsvFile Index_curves ForceIndex spring slope
FilePath = handles.file_path;
handles.Lc_Lp = [];
cla reset % reset figure
Index_curves = handles.Index_curves;
if handles.ForceIndex<length(Index_curves) && length(Index_curves)>1
    ForceIndex = handles.ForceIndex + 1;
else
    ForceIndex = handles.ForceIndex;
end

spring = handles.spring; % N/m
slope = handles.slope; % m/V

if isfield(handles, 'fit_curve')
    handles = rmfield(handles, 'fit_curve');
end
    
% AFM data
CsvFile = [FilePath 'Force' num2str(Index_curves(ForceIndex),'%03d') '.csv']; % [Extension Force] (nm, pN)
if exist(CsvFile,'file')
    data = dlmread(CsvFile,',',1,0);
    extension = data(:,1);
    force = data(:,2);
else
    data_retract_raw = handles.data_raw(:,2*(ForceIndex-1)+2:-1:2*(ForceIndex-1)+1); % [Extension Force]
    extension_raw = data_retract_raw(:,1); % unit: m
    force = data_retract_raw(:,2)*spring*slope*1e12; % unit: pN
    extension = extension_raw*1e9 - force/spring*1e-3; % unit: nm, consider the cantilever shape change
    
    % Correct the baseline and origin
    p = handles.p; % linear fit parameters
    p0 = handles.origin(:,:,ForceIndex);
    force = force - polyval(p, extension) - p0(1,2);
    extension = extension - p0(1,1);
end

% plot the retract curve
datacursormode on
baseX = min(extension):0.1:max(extension);
baseY = 0*baseX;
plot(baseX, baseY, '--r')
hold on
plot(extension, force)
xlabel('Extension (nm)');ylabel('Force (pN)');
title(['Force vs Extension (# ' num2str(Index_curves(ForceIndex)) ')'])
axis(handles.axisRangeZoomed);
hold off
handles.data_retract_cor = [extension force];
handles.data_retract = handles.data_retract_cor; % for smoothing data
handles.ForceIndex = ForceIndex;
guidata(hObject, handles);


function previous_curve_Callback(hObject, eventdata, handles)
global FilePath CsvFile Index_curves ForceIndex spring slope
FilePath = handles.file_path;
handles.Lc_Lp = [];
cla reset % reset figure
if handles.ForceIndex>1 && length(Index_curves)>1
    ForceIndex = handles.ForceIndex - 1;
else
    ForceIndex = handles.ForceIndex;
end

spring = handles.spring; % N/m
slope = handles.slope; % m/V

% AFM data
CsvFile = [FilePath 'Force' num2str(Index_curves(ForceIndex),'%03d') '.csv']; % [Extension Force] (nm, pN)
if exist(CsvFile,'file')
    data = dlmread(CsvFile,',',1,0);
    extension = data(:,1);
    force = data(:,2);
else
    data_retract_raw = handles.data_raw(:,2*(ForceIndex-1)+2:-1:2*(ForceIndex-1)+1); % [Extension Force]
    extension_raw = data_retract_raw(:,1); % unit: m
    force = data_retract_raw(:,2)*spring*slope*1e12; % unit: pN
    extension = extension_raw*1e9 - force/spring*1e-3; % unit: nm, consider the cantilever shape change
    % Correct the baseline and origin
    p = handles.p; % linear fit parameters
    p0 = handles.origin(:,:,ForceIndex);
    force = force - polyval(p, extension) - p0(1,2); % corrected
    extension = extension - p0(1,1);
end
% plot the retract curve
datacursormode on
baseX = min(extension):0.1:max(extension);
baseY = 0*baseX;
plot(baseX, baseY, '--r')
hold on
plot(extension, force)
xlabel('Extension (nm)');ylabel('Force (pN)');
title(['Force vs Extension (# ' num2str(Index_curves(ForceIndex)) ')'])
axis(handles.axisRangeZoomed);
hold off
handles.data_retract_cor = [extension force];
handles.data_retract = handles.data_retract_cor; % for smoothing data
handles.ForceIndex = ForceIndex;
guidata(hObject, handles);


function spring_value_Callback(hObject, eventdata, handles)
handles.spring = str2double(get(handles.spring_value,'string'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function spring_value_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function slope_value_Callback(hObject, eventdata, handles)
handles.slope = str2double(get(handles.slope_value,'string'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slope_value_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in baseline.
function baseline_Callback(hObject, eventdata, handles)
% hObject    handle to baseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global p ForceIndex
ForceIndex = handles.ForceIndex;
data = handles.data_retract_cor;
extension = data(:,1);
force = data(:,2);
x1 = handles.p1_x;
x2 = handles.p2_x;
index1 = find(extension >= x1 & extension<= x2);
extension1 = extension(index1);
force1 = force(index1);
p = polyfit(extension1,force1,1); % p = [k b]
force_cor = force - polyval(p, extension); % corrected

datacursormode on
baseX = min(extension):0.1:max(extension);
baseY = 0*baseX;
plot(baseX, baseY, '--r')
hold on
plot(extension,force_cor)
ylim([-150 max(force_cor)+50])
xlabel('Extension (nm)');ylabel('Force (pN)');
hold off

handles.data_retract_cor = [extension force_cor];
handles.p = p; % save linear fit parameters
guidata(hObject, handles);


% --- Executes on button press in baseline_pts.
function baseline_pts_Callback(hObject, eventdata, handles)
% hObject    handle to baseline_pts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dcm_obj = datacursormode(gcf);
info = getCursorInfo(dcm_obj);
handles.p1_x = min([info(1).Position(1) info(2).Position(1)]);
handles.p2_x = max([info(1).Position(1) info(2).Position(1)]);

guidata(hObject, handles);


% --- Executes on button press in origin.
function origin_Callback(hObject, eventdata, handles)
% hObject    handle to origin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global FilePath CsvFile data Index_curves ForceIndex;
FilePath = handles.file_path;
data = handles.data_retract_cor;
ForceIndex = handles.ForceIndex;
Index_curves = handles.Index_curves;
xRange = get(gca, 'xlim');
yRange = get(gca, 'ylim');

dcm_obj = datacursormode(gcf);
info = getCursorInfo(dcm_obj);
p0 = info.Position;

extension = data(:,1) - p0(1,1);
force = data(:,2) - p0(1,2);
ymin = -150;
extension = extension(force >= ymin); % filter the data
force = force(force >= ymin);
plot(extension, force)
xlabel('Extension (nm)');ylabel('Force (pN)');
title(['Force vs Extension (# ' num2str(Index_curves(ForceIndex)) ')'])
xLimits = xRange-p0(1,1); yLimits = yRange-p0(1,2);
xlim(xLimits)
ylim(yLimits)
hold on
baseX = min(extension):0.1:max(extension);
baseY = 0*baseX;
plot(baseX, baseY, '--g')
hold off

XSpan = xLimits(2)-xLimits(1); YSpan = yLimits(2)-yLimits(1);
set(handles.x_min,'string',num2str(xLimits(1),'%.0f'));
set(handles.x_max,'string',num2str(xLimits(2),'%.0f'));
set(handles.y_min,'string',num2str(yLimits(1),'%.0f'));
set(handles.y_max,'string',num2str(yLimits(2),'%.0f'));
set(handles.x_span,'string',num2str(XSpan,'%.0f'));
set(handles.y_span,'string',num2str(YSpan,'%.0f'));

% update data
handles.data_retract_cor = [extension force];
handles.origin(:,:,ForceIndex) = handles.origin(:,:,ForceIndex)+p0;
for i=1:length(Index_curves) % update the origin for other forces
    if any(handles.origin(:,:,i),2)==0 % check if all elements are 0
        handles.origin(:,:,i) = p0;
    end
end
handles.xLimits = xLimits;
handles.yLimits = yLimits;
handles.axisRangeZoomed = [xLimits yLimits];

% save data
CsvFile = [FilePath 'Force' num2str(Index_curves(ForceIndex),'%03d') '.csv']; % [Extension Force] (nm, pN)
%write header to file
header = 'Extension (nm), Force (pN)';
fid = fopen(CsvFile,'w');
fprintf(fid,'%s\n',header);
fclose(fid);
%write data to end of file
dlmwrite(CsvFile,handles.data_retract_cor,'-append', 'precision', '%.3f');

guidata(hObject, handles);


% --- Executes on button press in wlc_fit.
function wlc_fit_Callback(hObject, eventdata, handles)
% hObject    handle to wlc_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global kB T fit_max_force;
% if isfield(handles, 'fit_curve')
%     cla(handles.fit_curve)
% end

data =  handles.data_retract_cor;
extension = data(:,1);
force = data(:,2);

dcm_obj = datacursormode(gcf);
info = getCursorInfo(dcm_obj);
handles.p1_x = min([info(1).Position(1) info(2).Position(1)]);
handles.p2_x = max([info(1).Position(1) info(2).Position(1)]);

x1 = handles.p1_x;
x2 = handles.p2_x;

index2 = find(extension >= x1 & extension <= x2);
x = extension(index2);
y_force = force(index2);

% i = handles.PeakNum;
% define WLC model function (anonymous)
WLC_model = @(L, x)(1e12*kB*T/(L(2)*1e-9)*(1/4*(1-x/L(1)).^(-2)-1/4+x/L(1))); % L=[Lc Lp]
Lc_Lp_trial = [max(x)+0.1 0.1]; % Lc_Lp = [Lc0 Lp0] initial values
my_opts = statset('MaxIter', 100000); % number of max iteration time
[ L, R, J, CovB, MSE ] = nlinfit(x, y_force, WLC_model, Lc_Lp_trial, my_opts);

Lc = L(1); Lp = L(2);
Lc_Lp_temp = L; % temporarily save Lc_Lp

if ~isfield(handles, 'Lc_Lp')
    Lc_Lp = [];
else
    Lc_Lp = handles.Lc_Lp;
end

DisplayFormat = '%.2f'; % display format of values in GUI
set(handles.Lc_cur,'string',num2str(Lc_Lp_temp(1), DisplayFormat));
set(handles.Lp_cur,'string',num2str(Lc_Lp_temp(2), DisplayFormat));
if size(Lc_Lp,1)>=1
    Lc_Lp_pre = Lc_Lp(end,:);
    dLc_dLp = Lc_Lp_temp-Lc_Lp_pre;
    set(handles.Lc_pre,'string',num2str(Lc_Lp_pre(1), DisplayFormat));
    set(handles.Lp_pre,'string',num2str(Lc_Lp_pre(2), DisplayFormat));
    set(handles.dLc,'string',num2str(dLc_dLp(1), DisplayFormat));
else
    set(handles.Lc_pre,'string','NaN');
    set(handles.Lp_pre,'string','NaN');
    set(handles.dLc,'string','NaN');
end

% remove data cursor points on figure
delete(findall(gcf,'Type','hggroup'));
% plot the fitting results
x_fit = 0:0.02:(x2+50); % unit: nm
force_fit = 1e12*kB*T/(Lp*1e-9)*(1/4*(1-x_fit/Lc).^(-2)-1/4+x_fit/Lc);
% max_force = min(500,max(force)+50); % unit: pN
% max_force = min(500,max(y_force)+20); % unit: pN
fit_max_force = handles.fit_max_force; % unit: pN

index3 = find(force_fit >= fit_max_force, 1, 'first');
x_fit = x_fit(1:index3);
force_fit = force_fit(1:index3);
hold on
handles.fit_curve = plot(x_fit,force_fit,'-r','LineWidth',2);
hold off

ColumnRange = {'A1:B' 'C1:D' 'E1:F' 'G1:H' 'I1:J' 'K1:L' 'M1:N' 'O1:P' 'Q1:R' 'S1:T'};
handles.data_fit = [x_fit; force_fit]';
% save guioutput
warning('off', 'MATLAB:xlswrite:AddSheet');
% xlswrite(ExcelFile,handles.data_fit,'Fitted curves',[ColumnRange{i} num2str(size(handles.data_fit(:,1),1))]);

handles.PeakNum = handles.PeakNum+1;
handles.Lc_Lp_temp = Lc_Lp_temp;
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function Lc_pre_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function Lc_cur_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function dLc_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function Lp_pre_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function Lp_cur_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in update.
function update_Callback(hObject, eventdata, handles)
% hObject    handle to update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Index_curves ForceIndex spring slope
ForceIndex = handles.ForceIndex;
spring = handles.spring; % N/m
slope = handles.slope; % m/V

% AFM data
data_retract_raw = handles.data_raw(:,2*(ForceIndex-1)+2:-1:2*(ForceIndex-1)+1); % [Extension Force]
extension_raw = data_retract_raw(:,1); % unit: m
force_raw = data_retract_raw(:,2)*spring*slope*1e12; % unit: pN
extension_cor = extension_raw*1e9 - force_raw/spring*1e-3; % unit: nm, consider the cantilever shape change
% plot the retract curve
datacursormode on
plot(extension_cor, force_raw)
xlabel('Extension (nm)');ylabel('Force (pN)');
title(['Force vs Extension (# ' num2str(Index_curves(ForceIndex)) ')'])
axis(handles.axisRangeZoomed);

handles.data_retract_cor = [extension_cor force_raw];
handles.data_retract = handles.data_retract_cor; % for smoothing data
handles.ForceIndex = ForceIndex;
guidata(hObject, handles);


% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close all;

% --- Executes on button press in save_fit.
function save_fit_Callback(hObject, eventdata, handles)
% hObject    handle to save_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Index_curves ForceIndex ForceNum Lc_Lp;
Index_curves = handles.Index_curves;
ForceIndex = handles.ForceIndex;
ForceNum = Index_curves(ForceIndex);
Lc_Lp_unsort = handles.Lc_Lp;
Lc_Lp = sortrows(Lc_Lp_unsort,1); % sort by contour length
filename = 'WLC-param.csv'; % [Extension Force] (nm, pN)

if exist(filename,'file')
    fid = fopen(filename,'a+');
else
    %write header to file
    fid = fopen(filename,'w+');
    header = 'Force #, Lc (nm), Lp (pN), dLc (nm)';    
    fprintf(fid,'%s\r\n',header);
end

for ii = 1:size(handles.Lc_Lp,1)
    if ii > 1
        fprintf(fid,'%3d,%.6f,%.6f,%.6f\r\n',ForceNum,Lc_Lp(ii,1),Lc_Lp(ii,2),Lc_Lp(ii,1)-Lc_Lp(ii-1,1));
%         fprintf(fid, '%9.5f  %8.5f  %8.5f\r\n', handles.Lc_Lp(ii,1),handles.Lc_Lp(ii,2),handles.Lc_Lp(ii,1)-handles.Lc_Lp(ii-1,1));
    else
        fprintf(fid,'%3d,%.6f,%.6f,%.6f\r\n',ForceNum,Lc_Lp(ii,1),Lc_Lp(ii,2),NaN);
%         fprintf(fid, '%9.5f  %8.5f  %8.5f\r\n', handles.Lc_Lp(ii,1),handles.Lc_Lp(ii,2),0);
    end
end
fclose(fid);

% create a new figure for saving
fig = figure('visible','off');
% copy axes into the new figure
ax3 = copyobj(handles.axes1,fig);
set(ax3, 'units', 'normalized', 'position', [0.1 0.2 0.8 0.6]);

% saveas(fig,'Force vs Extension.png')
print(fig,'-depsc',[handle.fig_path 'Force vs Extension - ' num2str(ForceNum,'%03d') '.eps'])
close(fig);


function filename_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function filename_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function unit_spring_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function unit_spring_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function unit_slope_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function unit_slope_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function update_range_Callback(hObject, eventdata, handles)
xLimits(1,1) = str2double(get(handles.x_min,'String'));
xLimits(1,2) = str2double(get(handles.x_max,'String'));
yLimits(1,1) = str2double(get(handles.y_min,'String'));
yLimits(1,2) = str2double(get(handles.y_max,'String'));
xlim(xLimits)
ylim(yLimits)

handles.xLimits = xLimits;
handles.yLimits = yLimits;
guidata(hObject, handles);


function x_span_Callback(hObject, eventdata, handles)
xLimits = get(gca,'XLim');
XSpan = str2double(get(handles.x_span,'String'));
xLimits(1,2) = xLimits(1,1)+XSpan;
set(handles.x_max,'string',num2str(xLimits(2),'%.0f'));
xlim(xLimits)

handles.XSpan = XSpan;
handles.xLimits = xLimits;
guidata(hObject, handles);

function x_span_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function y_span_Callback(hObject, eventdata, handles)
yLimits = get(gca,'YLim');
YSpan = str2double(get(handles.y_span,'String'));
yLimits(1,2) = yLimits(1,1)+YSpan;
set(handles.y_max,'string',num2str(yLimits(2),'%.0f'));
ylim(yLimits)

handles.YSpan = YSpan;
handles.yLimits = yLimits;
guidata(hObject, handles);

function y_span_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function x_min_Callback(hObject, eventdata, handles)
xLimits = get(gca,'XLim');
xLimits(1,1) = str2double(get(handles.x_min,'String'));
xLimits(1,2) = xLimits(1,1)+handles.XSpan;
set(handles.x_max,'string',num2str(xLimits(2),'%.0f'));
xlim(xLimits)

handles.xLimits = xLimits;
guidata(hObject, handles);

function x_min_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function x_max_Callback(hObject, eventdata, handles)
xLimits = get(gca,'XLim');
xLimits(1,2) = str2double(get(handles.x_max,'String'));
xLimits(1,1) = xLimits(1,2)-handles.XSpan;
set(handles.x_min,'string',num2str(xLimits(1),'%.0f'));
xlim(xLimits)

handles.xLimits = xLimits;
guidata(hObject, handles);

function x_max_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function y_min_Callback(hObject, eventdata, handles)
yLimits = get(gca,'YLim');
yLimits(1,1) = str2double(get(handles.y_min,'String'));
yLimits(1,2) = yLimits(1,1)+handles.YSpan;
set(handles.y_max,'string',num2str(yLimits(2),'%.0f'));
ylim(yLimits)

handles.yLimits = yLimits;
guidata(hObject, handles);

function y_min_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function y_max_Callback(hObject, eventdata, handles)
yLimits = get(gca,'YLim');
yLimits(1,2) = str2double(get(handles.y_max,'String'));
yLimits(1,1) = yLimits(1,2)-handles.YSpan;
set(handles.y_min,'string',num2str(yLimits(1),'%.0f'));
ylim(yLimits)

handles.yLimits = yLimits;
guidata(hObject, handles);

function y_max_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function get_zoomed_range_Callback(hObject, eventdata, handles)
xLimits = get(gca,'XLim');
yLimits = get(gca,'YLim');

handles.axisRangeZoomed = [xLimits yLimits];
guidata(hObject, handles);


function draw_baseline_Callback(hObject, eventdata, handles)
global ForceIndex
ForceIndex = handles.ForceIndex;
data = handles.data_retract_cor;
extension = data(:,1);
force = data(:,2);
x1 = handles.p1_x;
x2 = handles.p2_x;
index1 = find(extension >= x1 & extension<= x2);
force1 = force(index1);
mean_force = mean(force1);

% datacursormode on
hold on
baseX = min(extension):0.1:max(extension);
baseY = mean_force*ones(1,length(baseX));
plot(baseX, baseY, '--g')
delete(findall(gcf,'Type','hggroup'));
hold off


function auto_scale_Callback(hObject, eventdata, handles)
axis tight


function fit_good_Callback(hObject, eventdata, handles)
if ~isfield(handles, 'Lc_Lp') || isempty(handles.Lc_Lp)
    handles.Lc_Lp = [handles.Lc_Lp_temp]; % assign the good fitting parameters
else
    if ~ismember(handles.Lc_Lp_temp,handles.Lc_Lp,'rows') % if Lc_Lp_temp not in Lc_Lp
        handles.Lc_Lp = [handles.Lc_Lp; handles.Lc_Lp_temp]; % append good fitting parameters
    end
end
guidata(hObject, handles);

function y_max_fit_Callback(hObject, eventdata, handles)
global ForceIndex fit_max_force;
ForceIndex = handles.ForceIndex;
data = handles.data_retract_cor;
extension = data(:,1);

fit_max_force = str2double(get(handles.y_max_fit,'String'));
hold on
xx = 0:0.1:max(extension);
yy = fit_max_force*ones(1,length(xx));
plot(xx,yy,'--r')
hold off

handles.fit_max_force = fit_max_force;
guidata(hObject, handles);

function y_max_fit_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
