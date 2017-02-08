function [ PeakInfo ] = FindPeaks( time, distance, force, min_height1, min_height2, min_peak_dist )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[pks_max, locs_max] = findpeaks( distance, 'minpeakheight',  min_height1, 'minpeakdistance', min_peak_dist);
[pks_min, locs_min] = findpeaks(-distance, 'minpeakheight', -min_height2, 'minpeakdistance', min_peak_dist);
max_dist_point = [locs_max, time(locs_max), distance(locs_max), force(locs_max), ones(length(locs_max),1)];
min_dist_point = [locs_min, time(locs_min), distance(locs_min), force(locs_min), zeros(length(locs_min),1)];
m_dist_point = [max_dist_point; min_dist_point];
[B, IX] = sort(m_dist_point(:, 1)); % according to time order
PossiblePeaks = m_dist_point(IX,:); % re-arrange the peak points
% take into account the two endpoints as peak points
PossiblePeaks = [1 time(1) distance(1) force(1) distance(1)>PossiblePeaks(1,3); ...
    PossiblePeaks; length(time) time(end) distance(end) force(end) distance(end)>PossiblePeaks(end,3)]; 
    
% Group consective max peaks or min peaks
peak_delete = [];
SameGroup = [1];
for ii = 2:length(PossiblePeaks)
    if PossiblePeaks(ii,5) == PossiblePeaks(ii-1,5)
        SameGroup = [SameGroup, ii]; %#ok<*AGROW>
    else
        % find the largest peaks
        if PossiblePeaks(SameGroup(1),5) == 1
            [C1,max_ind] = max(PossiblePeaks(SameGroup,3));
            SameGroup(max_ind) = [];
            peak_delete = [peak_delete SameGroup];
        else
            % find the smallest peaks
            [C2,min_ind] = min(PossiblePeaks(SameGroup,3));
            SameGroup(min_ind) = [];
            peak_delete = [peak_delete SameGroup];
        end
        SameGroup = [ii];
    end
    if ii == length(PossiblePeaks)
        if PossiblePeaks(SameGroup(1),5) == 1
            [C1,max_ind] = max(PossiblePeaks(SameGroup,3));
            SameGroup(max_ind) = [];
            peak_delete = [peak_delete SameGroup];
        else
            [C2,min_ind] = min(PossiblePeaks(SameGroup,3));
            SameGroup(min_ind) = [];
            peak_delete = [peak_delete SameGroup];
        end
    end
end

[TruePeaks,PS] = removerows(PossiblePeaks,'ind',peak_delete); % remove points that are not real peaks
PeakInfo = TruePeaks(:,2:5); % [time distance force group], group = 0/1

end