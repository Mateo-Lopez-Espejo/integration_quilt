% Top-level script for generating the quilts for the ferrets.
% The code was adapted from synth_quilt_naturalsound10.m in the
% ecog-quilting/code directory.

% 2018-02-25: Last edited, Sam NH
% 
% 2018-03: Modified to use whitened stimuli, Sam NH

% paths to external code repositories
addpath(genpath([root_directory '/Sound_Quilting_Toolbox_v1.2']));
addpath(genpath([root_directory '/general-analysis-code']));
addpath(genpath([root_directory '/general-audio-code']));
addpath(genpath([root_directory '/sbatch-code']));

% top-level directory for this project
project_directory = [root_directory '/scrambling-ferrets'];

%% Parameters

clear P;

% sampling rate
P.audio_sr = 44096;

% rms level
P.rms_level = 0.01;

% duration of the smallest segment in milliseconds and samples
P.win_sec = 1/64;
P.win_ms = 1000*P.win_sec;
P.win_smp = checkint(P.audio_sr*P.win_sec);

% segment sizes
P.seg_sizes_in_win = [1,2,4,8,16,32];

% number of random orderings to create
P.n_orders = 2;

% number of total stimuli (must match number in input directory)
% and duration of sound segment to take from each stimulus
P.n_stims = 20;
P.duration_per_stim_sec = 0.5;
total_duration = P.duration_per_stim_sec * P.n_stims;

% 2 = window, but not PSOLA, 3 = PSOLA
P.concatenation_method = 2;

% total duration of a single quilt
% needs to evenly divide the total duration
% which is currently 10 seconds (20 stim * 0.5 sec each)
P.quilt_duration_sec = 10;
n_quilts = checkint(total_duration / P.quilt_duration_sec);

% duration in which to expect an onset response
% scrambling is done separately for this window
% to make it possible to just analyze the last N-onset_response_region
% seconds of the stimulus
P.onset_response_region = 0.5;

% whether or not to generate practice stimuli
practice = false;

%% Directories / stimuli

% stimuli and directories
source_directory = [root_directory ...
    '/scrambling-ferrets/stimuli/naturalsound-v2-whitened'];
if practice; source_directory = [source_directory '-practice']; end
quilt_directory = [source_directory '-quilt-' num2str(P.duration_per_stim_sec) 'sec' ...
    '-catmethod' num2str(P.concatenation_method)];
if ~exist(quilt_directory, 'dir');
    mkdir(quilt_directory);
end
stims = mydir(source_directory);
assert(length(stims)==P.n_stims);

%% Main

for k = 1:length(P.seg_sizes_in_win) % loop through seg sizes
    
    % segment size in seconds and samples
    seg_smp = checkint(P.win_smp * P.seg_sizes_in_win(k));
    seg_sec = seg_smp / P.audio_sr;
    fprintf('Seg size %.3f sec\n', seg_sec);
    drawnow;
    
    % duration of the stim plus buffer time
    % we need an extra segment of buffer on either end for the purposes
    % of the overlap and add method
    duration_plus_buffer_sec = P.duration_per_stim_sec + seg_sec*2;
    
    % number of segments per stimulus excluding the buffer
    n_segs_per_stim = checkint(P.duration_per_stim_sec / seg_sec);
    
    % number of segments per quilt
    n_segs_per_quilt = checkint(P.quilt_duration_sec / seg_sec);
    
    % load the source material and determine the segment numbers
    % note this is done separately for every segment size because of the segment
    % dependent buffer
    count = 0;
    source_allstims = nan(duration_plus_buffer_sec * P.audio_sr, P.n_stims);
    segs_allstims = nan(n_segs_per_stim, P.n_stims);
    for j = 1:P.n_stims
                
        % read in waveform and format
        [y,sr] = audioread([source_directory '/' stims{j}]);
        
        % resample to desired audio rate
        y = resample(y, P.audio_sr, sr);
        clear sr;
        
        % select first duration_per_stim_sec seconds
        y = y(1:P.duration_per_stim_sec*P.audio_sr);
                
        % rms normalize
        y = P.rms_level * y / rms(y(:));
        
        % save relevant bit with zero-padded buffer
        source_allstims(:,j) = [zeros(seg_smp,1); y; zeros(seg_smp,1)];
        
        % add 1 to account for the buffer
        segs_not_in_buffer = (1:n_segs_per_stim) + 1;
        n_segs_per_stim_including_buffer = (n_segs_per_stim+2);
        segs_allstims(:,j) = ...
            segs_not_in_buffer + n_segs_per_stim_including_buffer*(j-1);
        
    end
    
    % time region in which to expect an onset response
    n_segs_in_onset_region = ceil(P.onset_response_region / seg_sec);
    
    overwrite = false;
    MAT_file = [quilt_directory '/randorders-seg-' num2str(round(seg_sec*1000)) 'ms.mat'];
    if ~exist(MAT_file, 'file') || overwrite
        
        % fix the random seed
        ResetRandStream2(k);
        
        % initialize and create first ordering
        rand_seg_order = nan(n_segs_per_quilt, n_quilts, P.n_orders);
        rand_seg_order(:,:,1) = reshape(...
            segs_allstims(randperm(numel(segs_allstims))),...
            [n_segs_per_quilt, n_quilts]);
        
        % create the next orders by randomizing within the stimulus
        for i = 2:P.n_orders
            for j = 1:n_quilts
                
                % separately shuffle onset and post onset
                while true
                    xi_onset = 1:n_segs_in_onset_region;
                    xi_post_onset = n_segs_in_onset_region + 1:n_segs_per_quilt;
                    xi_shuff = [Shuffle(xi_onset), Shuffle(xi_post_onset)];
                    
                    prev1 = nan(1, n_segs_per_quilt);
                    prev2 = nan(1, n_segs_per_quilt);
                    prev1(2:n_segs_per_quilt) = 1:n_segs_per_quilt-1;
                    prev2(xi_shuff(2:end)) = xi_shuff(1:end-1);
                    
                    yi = ~isnan(prev1) & ~isnan(prev2);
                    if all(prev1(yi) ~= prev2(yi))
                        fprintf('succeed\n')
                        break;
                    else
                        fprintf('fail\n')
                    end
                    clear yi;
                end
                
                rand_seg_order(:,j,i) = rand_seg_order(xi_shuff,j,1);
                
            end
        end
        
        save(MAT_file, 'rand_seg_order');
        
    else
        
        load(MAT_file, 'rand_seg_order');
    end
    
    % create the quilts and write them to file
    for i = 1:P.n_orders
        for j = 1:n_quilts
            
            % create quilt
            overlap_smp = P.win_smp;
            max_shift_smp = floor(P.win_smp/2);
            order = rand_seg_order(:,j,i);
            [quilt,shifts] = reorder_segs(...
                source_allstims(:), seg_smp, order, ...
                P.concatenation_method, overlap_smp, max_shift_smp);
                        
            %             % verify that things worked assuming no shift
            %             xi = 100;
            %             y = quilt((1:seg_smp) + (xi-1)*seg_smp);
            %             z = source_allstims((1:seg_smp) + (order(xi)-1)*seg_smp)';
            %             plot([y,z]);
                        
            % write to .wav file
            fname = [quilt_directory '/seg-' num2str(round(seg_sec*1000)) ...
                'ms-stim' num2str(j) '-order' num2str(i)];
            audiowrite_checkclipping([fname '.wav'], quilt, P.audio_sr);
            
            % write to mat file
            save([fname '.mat'], 'quilt', 'P', 'shifts', 'order', ...
                'n_segs_per_stim', 'n_segs_in_onset_region', 'n_segs_per_quilt');
                        
        end
    end
end


%% Check the orders

% for j = 34
% 
% k = 2;
% 
% j
% stims{j}
% % segment size in seconds and samples
% seg_smp = P.win_smp * P.seg_sizes_in_win(k);
% seg_sec = seg_smp / P.audio_sr;
% n_segs_per_stim = P.duration_per_stim_sec / seg_sec;
%     
% 
% [y,sr] = audioread([source_directory '/' stims{j}]);
% y = format_wav(y, sr, P);
% 
% % subdivide
% y_seg = reshape(y, [seg_smp, length(y)/seg_smp]);
% 
% % power 
% half1 = round(P.win_smp/2);
% half2 = P.win_smp - half1;
% win = [.5*(1-cos(pi*[0:half1-1]'/half1)); ones(seg_smp-P.win_smp,1); ...
%     .5*(cos(pi*[0:half2-1]'/half2)+1)];
% rms_pow = 10*log10(mean(bsxfun(@times, win, y_seg).^2,1));
% rms_pow = rms_pow - max(rms_pow);
% 
% 
% plot(rms_pow);
% hold on;
% plot([1, length(rms_pow)], [-30, -30], 'r--');
% sum(rms_pow > -30)*seg_sec
% end
% 
% % plot(y)
% %%
% 
% % remove breaths and pauses
% win_ms = 30;
% amp_thresh = -35;
% gap_len_thresh = 8;
% [y_dense, gap_contents] = excise_breaths_and_pauses(y,P.audio_sr,win_ms,amp_thresh,gap_len_thresh);
% % y_dense = 0.05 * y_dense / rms(y_dense);
% % length(gap_contents)/length(y)
% 
% plot(y)
% sound(y_dense, P.audio_sr)
% sound(y, P.audio_sr)