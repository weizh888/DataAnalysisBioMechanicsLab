clear variables
clc

global date_format startInd existExtensionCol
% Define Date format corresponding to the folder names.
date_format = 1; % Folder name: 2016-01-01
% date_format = 0; % Folder name:20160101

% For different naming method, such as '00 F vs LR_v200_unfold.txt' and 00 v200_unfold.txt
% startInd = 12; % '00 F vs LR_v200_unfold.txt'
startInd = 4; % '00 v200_unfold.txt'
speedRange = startInd:startInd+5;
existExtensionCol = 0; % If Extension column exists, remove it; else, don't remvove.

% Get a list of all files and folders in this folder.
folder_path = 'E:\00 Research\00 Projects\Project - GPIb-IX\Analysis\GPIb-IX WT_enzyme vs WM23';
listing = dir(folder_path);
listing(1:2) = [];
% Get a logical vector that tells which is a directory.
dirFlags = [listing.isdir];
% Extract only those that are directories.
subFolders = listing(dirFlags);
% Print folder names to command window.
for k = 1 : length(subFolders)
    % Extract experiment date from folder name.
    % Two different folder naming methods used.
    d = subFolders(k).name;
    switch date_format
        case 1 % Date: 2017-01-31
            date = [d(6:7) '/' d(9:10) '/' d(1:4)];
            % date = datestr([d(1:4) '/' d(6:7) '/' d(9:10)]);
        case 0 % Date: 20170131
            date = [d(5:6) '/' d(7:8) '/' d(1:4)];
            % date = datestr([d(1:4) '/' d(5:6) '/' d(7:8)]);
    end
    
    %% Get data files.
    txtListing = dir(fullfile([folder_path '/' subFolders(k).name], '*.txt'));
    file_list = {txtListing.name}'
    % Find data files, whose names contain 'unfold', 'refold' and
    % 'rupture'.
    expr = ('unfold|refold|rupture');
    ix=regexp(file_list,expr);
    ix=~cellfun('isempty',ix);
    file_list_new=file_list(ix)
    
    %% Check the existence of correction files.
    f_cor_path = [folder_path '/' subFolders(k).name '/Force_Cor_Factor.txt'];
    d_cor_path = [folder_path '/' subFolders(k).name '/Distance_Cor_Factor.txt'];
    if exist(f_cor_path,'file')==2
        fid0 = fopen(f_cor_path, 'r');
        data_f_cor = textscan(fid0, '%s %f', 'HeaderLines', 0, 'CollectOutput', 1);
        % Remove NaN rows.
        idx1 = find(isnan(data_f_cor{2}));
        data_f_cor{1}(idx1,:) = [];
        data_f_cor{2}(idx1,:) = [];
    else
        data_f_cor{1} = 'NA';
        data_f_cor{2} = 1;
    end
    if exist(d_cor_path,'file')==2
        fid1 = fopen(d_cor_path, 'r');
        data_d_cor = textscan(fid1, '%s %f %f', 'HeaderLines', 0, 'CollectOutput', 1);
        idx2 = find(isnan(data_d_cor{2}));
        data_d_cor{1}(idx2,:) = [];
        data_d_cor{2}(idx2,:) = [];
    else
        data_f_cor{1} = 'NA';
        data_f_cor{2} = [1, 0];
    end
    
    %% Process data files
    if length(file_list_new)>=1 % If data file exists.
        for j=1:length(file_list_new)
            data_path = [folder_path '/' subFolders(k).name '/' file_list_new{j}];
            fid = fopen(data_path, 'r');
            data = textscan(fid, '%f %f %f %f %f %f %s', 'HeaderLines', 1, 'CollectOutput', 1);
            % Get the pulling speed.
            switch file_list_new{j}(speedRange) % (12:15)
                case 'v100'
                    speed = '100';
                    typeRange = speedRange+5; % 17:22
                case 'v200'
                    speed = '200';
                    typeRange = speedRange+5;
                case 'v400'
                    speed = '400';
                    typeRange = speedRange+5;
                case 'v500'
                    speed = '500';
                    typeRange = speedRange+5;
                case 'v750'
                    speed = '750';
                    typeRange = speedRange+5;
                case 'v300'
                    speed = '300';
                    typeRange = speedRange+5;
                case 'v250'
                    speed = '250';
                    typeRange = speedRange+5;
                case 'v150'
                    speed = '150';
                    typeRange = speedRange+5;
                case 'v50_'
                    speed = '50';
                    typeRange = speedRange+4; % Move one letter forward.
            end
            
            % get the event type: unfold, refold, or rupture
            switch file_list_new{j}(typeRange)
                case 'unfold'
                    type = 'unfold';
                    if existExtensionCol==1
                        data{1}(:,4) = []; % Remove one column, i.e., the extension column.
                    end
                case 'refold'
                    type = 'refold';
                    if existExtensionCol==1
                        data{1}(:,4) = [];
                    end
                case 'ruptur'
                    type = 'rupture';
            end
            
            disp([d ': ' type ', ' num2str(speed) 'nm/s'])
            write_to_file(speed, type, date, data, data_f_cor, data_d_cor);
        end
    end
end
fclose all;