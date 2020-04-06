function [L, M, raster_separate_reps, t, S] = ...
    lag_corr_cross_segdur(recording_directory, recording_id, varargin)

% Primary analysis script for performing the lagged correlation analysis.

% 2018-03-26: Created, Sam NH
% 
% 2019-04-15: Last updated, Sam NH

%% Parameters

% directories to save results to
I.figure_directory = pwd;
I.analysis_directory = pwd;
I.raster_directory = pwd;

I.lag_win = [0, 1];
I.plot_win = [0, 0.5];
I.boundary = 'none';
I.trancorr = 'none';
I.tranweight = 'none';
I.divnorm = true;
I.tranweightnsegs = 'none';
I.tranweightdenom = 'none';
I.weightdenom = true;
I.distr = 'gamma';
I.plot_figure = true;

% spike threshold
I.spike_thresh = -4;

% amount by which to smooth each spike
I.fwhm_ms = 10;

% sampling rate of the raster
I.raster_sr = 500;

% individual units to plot
I.units = [];

% whether or not to enter debug mode
I.keyboard = false;

% whether or not to use single unit data
I.single_unit = false;

% whether or not to overwrite results
I.overwrite = false;

% save formated rasters, by request of Sam. MLE
I.save_rasters = false; 

I = parse_optInputs_keyvalue(varargin, I);

if I.keyboard
    keyboard;
end

%% Parse stimulus info

T = read_trial_info(recording_directory, recording_id);

%% Compute raster

stim_dur_sec = 10;
win = [0, stim_dur_sec + T.prestim_silence*2];
[raster_smooth, tbins] = raster(recording_directory, recording_id, win, ...
    'fwhm_ms', I.fwhm_ms, 'raster_sr', I.raster_sr, 'spike_thresh', I.spike_thresh, ...
    'single_unit', I.single_unit);

n_units = size(raster_smooth,3);
if isempty(I.units)
    I.units = 1:n_units;
end

% LBHB cell names, MLE 2019 06 11
if I.single_unit
    chnames = lbhb_cell_name(recording_directory, recording_id);
else
    chnames = cell(1, n_units);
    for i=1:n_units
        formatSpec = '%s-%02d-e';
        site_name = recording_id(1:7); %lbhb site code
        chnames{i} = sprintf(formatSpec, site_name, i);
    end
end
 
    
%% Separate out repetitions into a separate dimension

n_conditions = T.n_seg * T.n_orders * T.n_stims;
% check where some combination of order and segment length 
% have less repetitions and breake the code therefore....
if mod(T.n_trials, n_conditions) ~= 0
    % chekc allpossible conditions, i.e, combinations of segment duration and order
    conditions = struct;
    c_count = 1;
    for o = 1:length(T.unique_orders)
        for s = 1:length(T.unique_segs)
            for st = 1:length(T.unique_stims)
                ord = T.unique_orders(o);
                seg = T.unique_segs(s);
                stm = T.unique_stims(st);
                conditions(c_count).orders = ord;
                conditions(c_count).segs = seg;
                conditions(c_count).stims = stm;
                mask = and(and(T.segs == seg, T.orders == ord), T.stims == stm);
                conditions(c_count).count = sum(mask);
                c_count = c_count + 1;
            end
        end
    end
    % minimum number of repetitions given all different contidiont 
    min_reps = min([conditions.count]);
    
    % get the extra trials of some conditions to get all the the min number
    extra_trials = [];
    for e = 1:length(conditions)
        if conditions(e).count > min_reps
            condition = and(and(T.orders == conditions(e).orders,...
                                T.segs == conditions(e).segs),...
                                T.stims == conditions(e).stims);
            cond_trial = T.trial_order(condition);
            extra_trials = [extra_trials, cond_trial(min_reps+1:end)];
        end
    end
    
    % excludes extra trials and recout trial number
    good_trials = ~ismember(T.trial_order,extra_trials);
    T.n_trials = length(good_trials);
    
    [~, T.trial_order] = sort(T.trial_order(good_trials));
    T.segs = T.segs(good_trials);
    T.orders = T.orders(good_trials);
    T.stims = T.stims(good_trials);
    T.stim_labels = T.stim_labels(good_trials);
    raster_smooth =  raster_smooth(:, good_trials, :);
    
    % sorts the trialse so when reshaping dimensions make sense...
    tosort = [T.segs; T.orders; T.stims];
    [~, sortidx] = sortrows(tosort', [1,2,3]);
    
    T.trial_order = T.trial_order(sortidx);
    T.segs = T.segs(sortidx);
    T.orders = T.orders(sortidx);
    T.stims = T.stims(sortidx);
    T.stim_labels = T.stim_labels(sortidx);
    raster_smooth =  raster_smooth(:, sortidx, :);
    
    T.n_trials = length(T.trial_order);

end

n_reps_total = T.n_trials/n_conditions;

% separate out different trials, permute dims
% time/bins x (repetitions and conditions) x electrodes
% -> time/bins x repetitions x conditions x electrodes
% -> time/bins x conditions x electrodes x repetitions
n_bins = size(raster_smooth,1);
raster_separate_reps = reshape(raster_smooth(:,T.trial_order,:), [n_bins, n_reps_total, n_conditions, n_units]);
raster_separate_reps = permute(raster_separate_reps, [1 3 4 2]);

%% Create the stimulus structure needed for the lag analysis

% one field per condition (removes repetitions)
clear S;
f = {'stim_labels', 'segs', 'stims', 'orders'};
for i = 1:length(f)
    S.(f{i}) = T.(f{i})(T.trial_order);
    S.(f{i}) = reshape(S.(f{i}), [n_reps_total, n_conditions]);
    S.(f{i}) = S.(f{i})(1,:)';
end

% segment durations need to be exact
S.segs(S.segs==16) = 15.625;
S.segs(S.segs==31) = 31.25;
S.segs(S.segs==63) = 62.50;

% ordering of the segments
F = filesep;
S.stim_directory = [root_directory F 'scrambling-ferrets' F 'stimuli' F ...
                    'naturalsound-v2-whitened-quilt-0.5sec-catmethod2'];
for i = 1:length(S.stim_labels)
    S.segorder(i) = load([S.stim_directory F S.stim_labels{i} '.mat']);
end

% additional duration information
S.scramstim_dur = 10; % total duration of each scrambled stimulus
S.sourcestim_dur = 0.5; % duration of the source stimuli

%% saves raw rasters
if I.save_rasters
    % select only good cells (reliability > 0.1)
    D = raster_separate_reps(:,:,I.units,:);
    t = tbins' - T.prestim_silence;
    cell_id = chnames(I.units);
    filename = mkpdir(I.raster_directory);
    save(filename, 'D', 't', 'S', 'cell_id');
else
    %% Perform lag analysis
    disp([newline 'performig lag analysis'])
    t = tbins' - T.prestim_silence;
    L = lag_corr_cross_segdur_modular(...
        raster_separate_reps, t, S, ...
        'boundary', I.boundary, 'lag_win', I.lag_win, ...
        'trancorr', I.trancorr, 'tranweight', I.tranweight, ...
        'plot_win', I.plot_win, 'channels', I.units, ...
        'output_directory', I.analysis_directory, ...
        'figure_directory', I.figure_directory, ...
        'chnames', chnames, 'single_unit', I.single_unit,...
        'overwrite', I.overwrite, ...
        'plot_figure', I.plot_figure);

    %% Fit lag correlations with model
    disp([newline 'fittig lags with model'])
    M = modelfit_lagcorr_cross_segdur_modular(L, ...
        'overwrite', I.overwrite, 'divnorm', I.divnorm, ...
        'tranweightnsegs', I.tranweightnsegs, 'tranweightdenom', I.tranweightdenom, ...
        'weightdenom', I.weightdenom, 'distr', I.distr, ...
        'plot_figure', I.plot_figure);
end

end