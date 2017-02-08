% Wei Zhang 5/1/2014
clear all; clc; close all;

% source_dir = 'E:\Dropbox\Research\Optical Tweezers Matlab - Wei\program\data';
source_dir = pwd;
dest_dir = source_dir;
source_files = dir(fullfile(source_dir, '*.txt'));

%% save data to a structure
% initialize the structure
field = {'v25' 'v40' 'v50' 'v100' 'v150' 'v200' 'v250' 'v300' 'v350'};
for field_num = 1:length(field)
    data.unfold.(field{field_num}) = {};
    data.refold.(field{field_num}) = {};
    data.rupture.(field{field_num}) = {};
end

for i = 1:length(source_files)
    fid1 = fopen(fullfile(source_dir, source_files(i).name), 'r'); % (Loading-rate [N/s], Rupture force [N], Extension [m], Time [s])
    if isempty(strfind(source_files(i).name,'rupture')) % rupture file only has 6 columns, while unfold file has 7 columns
        alldata = textscan(fid1, '%f  %f  %f  %f  %f  %f  %s' , 'headerlines', 1);
        alldata = [num2cell([alldata{1:6}]) alldata{7}];
    else
        alldata = textscan(fid1, '%f  %f  %f  %f  %f  %s' , 'headerlines', 1);
        alldata = [num2cell([alldata{1:5}]) alldata{6}];
    end
    C = strsplit(source_files(i).name,'_');
    alldata(cellfun(@isempty,alldata)) = C(2);
    
    % unfold data
    if ~isempty(strfind(source_files(i).name,'unfold'))
        for field_num = 1:length(field)
            if ~isempty(strfind(source_files(i).name, field{field_num}))
                data.unfold.(field{field_num}) = [data.unfold.(field{field_num}); alldata];
            end
        end
        % refold data
    elseif ~isempty(strfind(source_files(i).name,'refold'))
        for field_num = 1:length(field)
            if ~isempty(strfind(source_files(i).name, field{field_num}))
                data.refold.(field{field_num}) = [data.refold.(field{field_num}); alldata];
            end
        end
        % rupture data
    elseif ~isempty(strfind(source_files(i).name,'rupture'))
        for field_num = 1:length(field)
            if ~isempty(strfind(source_files(i).name, field{field_num}))
                data.rupture.(field{field_num}) = [data.rupture.(field{field_num}); alldata];
            end
        end
    end
end
fclose('all');
%% write data for existed speeds
for field_num = 1:length(field)
    if ~isempty(data.unfold.(field{field_num}))
        xlswrite(fullfile(dest_dir,'Force vs LR (unfold).xlsx'), data.unfold.(field{field_num}), [field{field_num} 'nm_s']);
    end
    if ~isempty(data.refold.(field{field_num}))
        xlswrite(fullfile(dest_dir,'Force vs LR (refold).xlsx'), data.refold.(field{field_num}), [field{field_num} 'nm_s']);
    end
    if ~isempty(data.rupture.(field{field_num}))
        xlswrite(fullfile(dest_dir,'Force vs LR (rupture).xlsx'), data.rupture.(field{field_num}), [field{field_num} 'nm_s']);
    end
end

%% delete default sheet 1, sheet 2, sheet 3
excelFileName = {'Force vs LR (unfold).xlsx' 'Force vs LR (refold).xlsx' 'Force vs LR (rupture).xlsx'};
excelFilePath = dest_dir;
% excelFilePath = pwd; % Current working directory.
sheetName = 'Sheet'; % EN: Sheet, DE: Tabelle, etc. (Lang. dependent)
for FileNum = 1:length(excelFileName)
    if exist([excelFilePath '\' excelFileName{FileNum}],'file')
        % Open Excel file.
        objExcel = actxserver('Excel.Application');
        objExcel.Workbooks.Open(fullfile(excelFilePath, excelFileName{FileNum})); % Full path is necessary!
        % Delete sheets.
        try
            % Throws an error if the sheets do not exist.
            objExcel.ActiveWorkbook.Worksheets.Item([sheetName '1']).Delete;
            objExcel.ActiveWorkbook.Worksheets.Item([sheetName '2']).Delete;
            objExcel.ActiveWorkbook.Worksheets.Item([sheetName '3']).Delete;
        catch
            ; % Do nothing.
        end
        % Save, close and clean up.
        objExcel.ActiveWorkbook.Save;
        objExcel.ActiveWorkbook.Close;
        objExcel.Quit;
        objExcel.delete;
    end
end