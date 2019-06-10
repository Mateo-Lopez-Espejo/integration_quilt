% Creates the stimulus orders for the Ferret scrambling experiment. The
% corresponding files are here (two directories contain the same stims):
% 
% /Users/svnh2/Dropbox (MIT)/mindhive/ecog-quilting/stimuli/naturalsound20-quilt-0.5sec-catmethod2
% 
% The stimuli are designed so that repetitions without reordering are spaced
% equally close as repetitions with reordering. 
% 
% 2018-02-25: Created, Sam NH

addpath(genpath([root_directory '/general-analysis-code']));

%% Generates sequences as though there were 8 runs

% segment durations
seg_durs = {'16', '31', '63', '125', '250', '500'};
n_seg_durs = length(seg_durs);

% total number of runs
n_blocks = 4;

%% stim pool

n_stims = 1;
n_orders = 2;
stim_pool = cell(n_seg_durs, n_stims, n_orders);
for i = 1:n_seg_durs
    for j = 1:n_stims
        for k = 1:n_orders
            
            % file name
            fname = ['seg-' seg_durs{i} 'ms-stim' num2str(j) ...
                '-order' num2str(k)];
            
            % assign
            stim_pool{i, j, k} = fname;
            clear fname index;
        
        end
    end
end

% unwrap
stim_pool = stim_pool(:);

%% randomize order within a pool per run

for r = 1:30
    
    ResetRandStream2(r);
    stim_order = cell(length(stim_pool), n_blocks);
    for q = 1:n_blocks
        stim_order(:,q) = stim_pool(randperm(length(stim_pool)));
    end
    
    stim_order = stim_order(:)
    
    output_directory = [root_directory '/scrambling-ferrets/stimuli/stim-orders20'];
    if ~exist(output_directory, 'dir'); mkdir(output_directory); end
    fname = [output_directory '/r' num2str(r) '.mat'];
    if ~exist(fname, 'file')
        save(fileplus(fname), 'stim_order');
    end
    
end
