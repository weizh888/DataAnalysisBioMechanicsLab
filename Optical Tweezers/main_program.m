% Wei Zhang (wez311@lehigh.edu)
% Lehigh University

% clearvars -except filename data_files
clear all;
close all; clc;

global filename file_path threshold_force
global UsefullData ExpDate

filename = input('Please input the data file to analyze: \n', 's');
file_path = ['data/' ExpDate '/' filename '/']; % set file path
ComFile = ['data/' filename 'COM.txt'];
if exist(ComFile, 'file')    
    ExpDate = input('Please input the date of experiment: \n', 's');
    if ~exist([pwd '/data/' ExpDate '/' filename],'dir') % create a new subfolder, if it doesn't exist
        mkdir([pwd '/data/' ExpDate], filename);
    end

    step1_load_data;
    
    if UsefullData == 1 % UsefullData is generated in step1_load_data
        step1_correct_dist_GUI;
        uiwait(step1_correct_dist_GUI);
        step2_find_peaks;
    end
    %%
    step3_split_curves;
    %
    threshold_force = 5; % unit: pN
    step4_select_curves;
    %
    step5_GUI;
    % Step_5_Unfolds_Rupture_new
    % Step_5_DNA_OS_Time
    % Step_5_High_Force_Hold_Time
    % Step5_GUI(filename,SpeedValue, style, threshold_force, Index_curves, si, PeakInfo)
    % Step_5_Plot_All_Curves
    % Step_5_Plot_Single_Curves
    % Step_5_Unfolds_new
    % Step_5_Rupture_new
    % Step_5_Rupture_No_Origin;
    % Step_5_Rupture;
    % Step_5_Unfold_No_Origin;
    % Step_5_Unfold_Refold_No_Origin
    % Step_5_Unfold_Rupture;
    % Step_5_WLC_fit;
    % % zero_dist = -1030; % unit of nm
    % zero_dist = input('Please input the indent of distance (nm):\n');
    % Step_5_extension_corrected;
else
    disp('The file doesn''t exist. Please input another file.')
    return
end