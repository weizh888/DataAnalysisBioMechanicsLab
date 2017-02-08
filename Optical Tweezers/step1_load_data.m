% Wei Zhang (wez311@lehigh.edu)
% Lehigh University
% Note add the name of loaded files, before starting the program

close all; clc; % clear all;
clearvars -except filename
global filename file_path raw_file SpeedValue UsefullData ExpDate
%% read all files
if ~exist('filename','var') || isempty(filename)
    filename = input('Please input the data file to analyze: \n', 's');
end
% if ~exist('SpeedValue','var') || isempty(SpeedValue)
fid = fopen([pwd '/data/' filename 'COM.txt']);
info = textscan(fid, '%s', 2, 'headerlines', 76);
SpeedValue = str2double(info{1}{2});
%     SpeedValue = input('Please input the pulling speed (nm/s): \n');
fclose(fid);
% end

if ~exist([pwd '/data/' ExpDate '/' filename],'dir') % create a new subfolder, if it doesn't exist
    mkdir([pwd '/data/' ExpDate], filename);
end

file_path = ['data/' ExpDate '/' filename '/']; % set file path
raw_file = [file_path 'All_Time_Force_TrapDist_raw.txt']; % [time distance force]
%%
FileList = {'A' 'B' 'C' 'D' 'E' 'F'};
FileNum = 1;
raw_data = [];
while FileNum <= 6
    raw_data_file = ['data/' filename FileList{FileNum} '.txt'];
    if exist(raw_data_file,'file')==2
        % read raw data file
        fid0 = fopen(raw_data_file,'r');
        temp_raw_data = [];
        if FileNum == 1
            % 1st file (file A), skip 2 lines at the beginning
            temp_raw_data = textscan(fid0, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'HeaderLines', 2, 'CollectOutput', 1, 'commentStyle', '#');
        else
            % other file, skip 1 line at the beginning
            temp_raw_data = textscan(fid0, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'HeaderLines', 1, 'CollectOutput', 1, 'commentStyle', '#');
        end
        fclose(fid0);
        raw_data = [raw_data; cell2mat(temp_raw_data)];
        % raw_data = dlmread(raw_data_file);
        FileNum = FileNum + 1;
    else
        % filename doesn't contain 'A','B','C'..
        raw_data_file = ['data/' filename '.txt'];
        if exist(raw_data_file,'file')==2
            % read raw data file
            fid0 = fopen(raw_data_file,'r');
            temp_raw_data = textscan(fid0, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'HeaderLines', 2, 'CollectOutput', 1, 'commentStyle', '#');
            fclose(fid0);
            raw_data = [raw_data; cell2mat(temp_raw_data)];
            break;
        else
            break;
        end
    end
end
raw_data = raw_data(1:end-1,:); % remove the last row, since the data might not be completely saved

%% pre-process raw data
% Define the Column Data that are detected
CycleCount = raw_data(:,1); % Sequential cycle count from main controller PIC (4000 Hz); Just divide CycleCount with 1000 to get the time in second.
Time = CycleCount./1000;
A_Psdy = raw_data(:,2);
A_LeverX = raw_data(:,3);
A_LeverY = raw_data(:,4);
B_Psdy = raw_data(:,5);
B_LeverX = raw_data(:,6);
B_LeverY = raw_data(:,7);
Y_force = raw_data(:,8);  % unit of pN
Tension = raw_data(:,9);  % unit of pN
A_dist_Y = raw_data(:,10);% Calculated A trap displacement in X (nm)
B_dist_Y = raw_data(:,11);% Calculated B trap displacement in X (nm)

A_dist_Y = smooth(A_dist_Y,32,'moving');
B_dist_Y = smooth(B_dist_Y,32,'moving');

% Y_dist = (A_dist_Y+B_dist_Y)./2; % unit of nm
% %% trap B distance is abnormal
DistCorFactor = 1; % default value
Y_dist = (A_dist_Y + B_dist_Y)./2.0; % unit of nm

A_CycleCount = raw_data(:,12);  % Cycle count from A trap PIC (16 bit, 4000 Hz)
B_CycleCount = raw_data(:,13);  % Cycle count from B trap PIC (16 bit, 4000 Hz)
time0 = raw_data(:,14);  % Calculated time from communications PIC cycle count (s)
Status = raw_data(:,15);  % Indicates state of pulling protocols (i.e., waiting, hopping, pulling, stretching)

%% plots of Dist_Y vs. Time, Y_Forve vs. Time, and Y_force vs.Dist
figure;
scnsize = get(0,'ScreenSize');
set(gcf, 'OuterPosition', [0 scnsize(4)*1/10 scnsize(3) scnsize(4)*9/10]);
h(1) = subplot(2,2,1);
plot(Time,Y_dist,'r');
ylabel('Distance (nm)');
xlim([min(Time) max(Time)]);ylim([min(Y_dist) max(Y_dist)]);
set(gca,'YGrid','on')

h(2) = subplot(2,2,3);
plot(Time,-Y_force,'g');
xlabel('Time (s)');
ylabel('Force (pN)');
xlim([min(Time) max(Time)]);

h(3) = subplot(2,2,[2 4]);
plot(Y_dist,-Y_force);
title('Force vs Distance');
xlabel('Distance (nm)');
ylabel('Force (pN)');
xlim([min(Y_dist) max(Y_dist)]);

set(h, 'box', 'off');
set(gca, 'LooseInset', get(gca, 'TightInset')); % set tight margin

%% save raw data
fid2 = fopen(raw_file,'w');
fprintf(fid2,'%8.3f %6.2f %9.4f %9.4f\r\n',[Time'; -Y_force'; A_dist_Y'; B_dist_Y';]);
fclose(fid2);

%% Construct a questdlg with three options
choice = questdlg('Does this data file contain any good curve(s)?', 'Good Data?', 'Yes','No','Close','Yes');
% Handle response
switch choice
    case 'No'
        UsefullData = 0;
        % removes the folder and its contents, if no good curves
        rmdir(['data/' ExpDate '/' filename ],'s')
        FileNum = 1;
        while FileNum <= 6
            fileToDelete = ['data/' filename FileList{FileNum} '.txt'];
            if exist(fileToDelete, 'file')
                delete(fileToDelete)
            end
            FileNum = FileNum+1;
        end
        delete(['data/' filename 'COM.txt'])
        close all; clc;
        break;
    case 'Close'
        close all; break;
    case 'Yes'
        UsefullData = 1;
end
fclose('all');