% Wei Zhang (wez311@lehigh.edu)
% Lehigh University

close all; clc;
clearvars -except filename data_files SpeedValue Min_Height1 Min_Height2 minpeak_distance style threshold_force PeakInfo ExpDate
global filename file_path ExpDate

if ~exist('filename','var') || isempty(filename)
    filename = input('Please input the data file to analyze: \n', 's');
end

% Define style value: style=1, then 1st curve is Approach; style=0, then 1st curve is Retract
if PeakInfo(1,2) > PeakInfo(2,2)
    style = 1;
else
    style = 0;
end

file_path = ['data/' ExpDate '/' filename '/']; % set file path
new_file = [file_path 'All_Time_Dist_Force.txt']; % [time distance force]
peak_file = [file_path 'All_Time_Dist_Force_peaks.txt'];  % [time distacne force]
single_cycle_file = [file_path 'Dist_Force_Time ']; % one approaching/retracting process

% Read the Time_Distance_Force file
data_peaks = dlmread(peak_file);
new_data = dlmread(new_file);
time = new_data(:,1);
distance = new_data(:,2);
force = new_data(:,3);

% threshold_force = 5; % define a max force to decide whether to show the force-extension curve
Index_curves = [];
TotalCycleNum = floor((length(data_peaks)+1)/2);
for kk = 1:TotalCycleNum
    str1 = [single_cycle_file num2str(kk) 'A.txt'];
    str2 = [single_cycle_file num2str(kk) 'R.txt'];
    % only analyze the data when retract curves exist
    if exist(str2, 'file') == 2
        data_retract = dlmread(str2);
        if exist(str1, 'file') == 2 % if approach file exists
            data_approach = dlmread(str1);
        end
        if max(data_retract(:,2) >= threshold_force) % decide whether something was pulled
            Index_curves = [Index_curves; kk]; % Add useful curve index
        else
            if exist(str1, 'file')==2
                delete(str1);
            end
            if exist(str2, 'file')==2
                delete(str2);
            end
        end
    end
end

fid = fopen('InfoForStep5.txt','w');
fprintf(fid, 'Date:\n%s\n', ExpDate);
fprintf(fid, 'Filename:\n%s\n', filename);
fprintf(fid, 'Pulling Speed (nm/s):\n%.2f\n', SpeedValue);
fprintf(fid, 'Seclected Curves:\n');
fprintf(fid, '%d ', Index_curves);
fclose(fid);