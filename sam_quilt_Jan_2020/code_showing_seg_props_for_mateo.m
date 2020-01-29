% load S structure sent by Mateo
%load([project_directory '/data/lag_corr_rasters_from_mateo_2019-06-26/AMT026a.mat'], 'S');
load([project_directory '/data/lag_corr_rasters_from_mateo_2019-06-26/AMT026a.mat'], 'S');

% number of scrambled stimuli, should be 12
n_scram_stim = length(S.stim_labels);

% duration of each scrambled stimulus is 10 seconds
S.scramstim_dur = 10;

% duration of the original source stimuli used 
% from which the scrambled stimuli were created
S.sourcestim_dur = 0.5;

% total number of original source stimuli
% there were 20 500 ms source stimuli
n_sourcestim = checkint(S.scramstim_dur/S.sourcestim_dur);

% now let's loop through all 12 stimuli
for stim_index = 1:n_scram_stim
    
    % segment duration for this scrambled stimulus
    seg_dur = S.segs(stim_index)/1000;
    
    % number of segments in each source stimulus
    % i.e. 500 ms / segment duration
    n_segs_per_sourcestim = S.sourcestim_dur / seg_dur;
    
    % this is the key
    % seg_order_indices is a matrix indicate the position and source stimulus for each segment
    seg_order_indices = reshape(sort(S.segorder(stim_index).order), n_segs_per_sourcestim, n_sourcestim);
    
    % now let's use the above matrix to specify the source and onset time
    % for each segment
    S.segorder(stim_index).source_index = nan(size(S.segorder(stim_index).order));
    S.segorder(stim_index).source_onset_time = nan(size(S.segorder(stim_index).order));
    for j = 1:length(S.segorder(stim_index).order) % looping over all segments
        
        % find the segment in the seg_order_indices matrix
        xi = find(S.segorder(stim_index).order(j) == seg_order_indices);
        
        % figure out its row and column position in the matrix
        [row_pos, col_pos] = ind2sub(size(seg_order_indices), xi);
        
        % convert to the desired onset time and source index
        % the column position is the source index
        S.segorder(stim_index).source_index(j) = col_pos;
        S.segorder(stim_index).source_onset_time(j) = (row_pos-1)*seg_dur;
        
    end
end

%% Now let's check the above is right

% let's pick the 20th segment from a 31 ms scrambled stimulus
directory_with_scrambled_stimuli = '/Users/svnh2/Dropbox (MIT)/mindhive/scrambling-ferrets/stimuli/naturalsound-v2-whitened-quilt-0.5sec-catmethod2';
stim_index = 3;
seg_index = 20;
seg_dur = S.segs(stim_index)/1000;
[wav, sr] = audioread([directory_with_scrambled_stimuli '/' S.stim_labels{stim_index} '.wav']);
inds = (1:round(seg_dur*sr)) + round((seg_index-1)*seg_dur*sr);
seg_from_scrambled_stim = wav(inds);

% now let's try to find the corresponding segment in the source stimuli
% note mydir is a personal function
directory_with_source_stimuli = '/Users/svnh2/Dropbox (MIT)/mindhive/scrambling-ferrets/stimuli/naturalsound-v2-whitened';
S.sources = mydir(directory_with_source_stimuli, '.wav');
source = S.sources{S.segorder(stim_index).source_index(seg_index)};
source_onset_time = S.segorder(stim_index).source_onset_time(seg_index);
[wav, sr] = audioread([directory_with_source_stimuli '/' source]);
inds = (1:round(seg_dur*sr)) + round(source_onset_time*sr);
seg_from_source_stim = wav(inds);

% now let's plot both to confirm they are the same
figure;
plot([seg_from_source_stim, seg_from_scrambled_stim])

