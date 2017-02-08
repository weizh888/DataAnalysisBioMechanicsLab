function [] = write_to_file(speed, type, date, data, data_f_cor, data_d_cor)
% speed: 100, 200, 400
% type: unfold, refold, rupture
% date: string, like '2017-01-01'
% data: txt files
filename = ['summary - ' num2str(speed) '_' type '.csv'];
% check the header
if exist(filename,'file')==2 % if file exists
    fid = fopen(filename, 'a+'); % append new results
else
    fid = fopen(filename, 'w+');
    switch type
        case 'unfold'
            fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\r\n',...
                'LR(N/s)', 'Unfold_F(N)', 'Unfold_D(m)', 'Time(s)',...
                'Curve#', 'FileName', 'Date', 'F_cor', 'Std (F_cor)', 'Ext_cor', 'Std (Ext_cor)');
        case 'refold'
            fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\r\n',...
                'LR(N/s)', 'Refold_F(N)', 'Refold_D(m)', 'Time(s)',...
                'Curve#', 'FileName', 'Date', 'F_cor', 'Std (F_cor)', 'Ext_cor', 'Std (Ext_cor)');
        case 'rupture'
            fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\r\n',...
                'LR(N/s)', 'Rupture(N)', 'Rupture_D(m)', 'Time(s)',...
                'Curve#', 'FileName', 'Date', 'F_cor', 'Std (F_cor)', 'Ext_cor', 'Std (Ext_cor)');
    end
end

if strcmp(type,'refold')==1
    data{1} = abs(data{1});
end

for k=1:size(data{1},1)
    % Get correction factors
    % data{2}{k} is the filename
    idx1 = find(ismember(data_f_cor{1}, data{2}{k}));

    if isempty(idx1)
        f_cor = mean(data_f_cor{2}); % Use mean value to replace it.
    elseif length(idx1)>=2
        f_cor = data_f_cor{2}(idx1(end));
    else
        f_cor = data_f_cor{2}(idx1);
    end
    % f_cor_std is not calculated using all force correction factors
    f_cor_std = std(data_f_cor{2});
    
    idx2 = find(ismember(data_d_cor{1}, data{2}{k}));
    if isempty(idx2)
        d_cor = mean(data_d_cor{2});
        d_cor_std = std(data_d_cor{2});
    elseif length(idx2)>=2
        d_cor = data_d_cor{2}(idx2(end),1);
        d_cor_std = data_d_cor{2}(idx2(1),2);
    else
        d_cor = data_d_cor{2}(idx2,1);
        d_cor_std = data_d_cor{2}(idx2,2);
    end
    
    fprintf(fid,'%10.4e, %10.4e, %10.4e, %.3f, %3d, %s, %s, %.4f, %.4f, %.4f, %.4f\r\n',...
        data{1}(k,1),data{1}(k,2),data{1}(k,3),data{1}(k,4),data{1}(k,5),data{2}{k},...
        date,f_cor,f_cor_std,d_cor,d_cor_std);
end

fclose all;
end