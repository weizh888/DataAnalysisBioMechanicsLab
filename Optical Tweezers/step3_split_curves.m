% Wei Zhang (wez311@lehigh.edu)
% Lehigh University

close all; clc;
clearvars -except filename data_files SpeedValue Min_Height1 Min_Height2 minpeak_distance...
    PeakInfo threshold_force ExpDate
global filename file_path ExpDate

if ~exist('filename','var') || isempty(filename)
    filename = input('Please input the data file to analyze: \n', 's');
end

% style=1: 1st curve is Approach; style=0: 1st curve is Retract.
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

for k = 1:length(data_peaks(:,1))-1 % save each curve [distance(nm) force(pN) time(s)]
    Time_k = [];
    Dist_Y_k = [];
    Y_force_k = [];
    j = 1;
    for i = 1:length(new_data) % split curves
        if time(i) >= data_peaks(k,1) && time(i) < data_peaks(k+1,1)
            Time_k(j) = time(i);
            Dist_Y_k(j) = distance(i);
            Y_force_k(j) = force(i);
            j = j+1;
        end
    end
    % first cuve is approaching curve
    if style == 1
        if mod(k, 2) == 1 % odd #
            str = [single_cycle_file num2str(floor((k+1)/2)) 'A.txt']; % save each curve (approaching)
        else
            str = [single_cycle_file num2str(floor((k+2)/2)) 'R.txt']; % save each curve (retracting)
        end
    else % first curve is retracting curve
        if mod(k, 2) == 1 % odd #
            str = [single_cycle_file num2str(floor((k+1)/2)) 'R.txt']; % save each curve (approaching)
        else
            str = [single_cycle_file num2str(floor((k+1)/2)) 'A.txt']; % save each curve (retracting)
        end
    end
    % Do not save the first approach file if retract not exist.
    if style~=1 || (style==1 && k~=1)
        fid = fopen(str,'w');
        DistForceTime_k = [Dist_Y_k', Y_force_k',Time_k']';
        fprintf(fid,'%9.3f %6.2f %8.3f\r\n',DistForceTime_k);
        fclose(fid);
    end
end