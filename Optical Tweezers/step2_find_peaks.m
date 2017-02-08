% Wei Zhang (wez311@lehigh.edu)
% Lehigh University

close all; clc;
clearvars -except filename SpeedValue min_height1 min_height2 min_peak_dist flag ExpDate
global filename file_path flag ExpDate

if ~exist('filename','var') || isempty(filename)
    filename = input('Please input the data file to analyze: \n', 's');
end

file_path = ['data/' ExpDate '/' filename '/']; % set file path
new_file = [file_path 'All_Time_Dist_Force.txt']; % [time distance force]
peak_file = [file_path 'All_Time_Dist_Force_peaks.txt'];  % [time distacne force]

new_data = dlmread(new_file); % [time dist force]
time = new_data(:, 1);
distance = new_data(:, 2);
force = new_data(:, 3);

flag = 0; % flag indicates whether to select peak threshold
min_height1 = mean(distance)+1;
min_height2 = mean(distance)-1;
min_peak_dist = 20;
while flag == 0
    [ PeakInfo ] = FindPeaks( time, distance, force, min_height1, min_height2, min_peak_dist );
    PeakMax = PeakInfo(PeakInfo(:,4) == 1,:);
    PeakMin = PeakInfo(PeakInfo(:,4) == 0,:);
    
    % plot of Dist vs Time
    figure;
    scnsize = get(0,'ScreenSize');
    set(gcf, 'OuterPosition', [0 scnsize(4)*1/10 scnsize(3) scnsize(4)*9/10]);
    plot(time,distance,'b',PeakMax(:, 1),PeakMax(:, 2),'or',PeakMin(:, 1),PeakMin(:, 2),'om','MarkerSize',6);
    title('Distance vs Time')
    xlabel('Time (s)');
    ylabel('Distance (nm)');
    xlim([min(new_data(:, 1)) max(new_data(:, 1))]);
    set(gca,'YGrid','on')
    
    reply = questdlg('Are all peaks seclected?', 'Curve Splitting', 'Yes', 'No', 'Close', 'Yes');
    switch reply
        case 'Yes'
            flag = 1; % if all peaks are selected, then jump out to next step
        case 'No'
            flag = 0;
            % Choose proper values of Min_Height1, Min_Height2, and minpeak_distance
            min_height1 = input('the min height of high peaks:\n');
            min_height2 = input('the max height of low peaks:\n');
        case 'Close'
            close all;
            break;
    end
end

% Save Peaks Data [time distance force]
fid1 = fopen(peak_file, 'w');
fprintf(fid1, '%8.3f %9.3f %6.2f\r\n', PeakInfo(:,1:3)');
fclose(fid1);