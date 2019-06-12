function r = reliability(recording_directory, recording_id, varargin)

% Measures test-retest reliability for all electrodes in a given recording.

% 2018-03-18: Last edited, Sam NH

%% Parameters

I.figure_directory = '';

% spike threshold
I.spike_thresh = -4;

% amount by which to smooth each spike
I.fwhm_ms = 10;

% sampling rate of the raster
I.raster_sr = 500;

% number of reps to use
I.n_reps = NaN;

% individual electrodes to plot
I.individ_units = [];

% whether or not to enter debug mode
I.keyboard = false;

% whether or not to use single unit data
I.single_unit = false;

% whether to plot responses
I.plot = true;

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

%% Plot test retest reliability

n_conditions = T.n_seg * T.n_orders * T.n_stims;
n_reps_total = T.n_trials/n_conditions;

% use all reps
if isnan(I.n_reps)
    I.n_reps = n_reps_total;
end

% separate out different trials
% time/bins x repetitions x conditions x electrodes
n_bins = size(raster_smooth,1);
n_electrodes = size(raster_smooth,3);
raster_separate_reps = reshape(raster_smooth(:,T.trial_order,:), [n_bins, n_reps_total, n_conditions, n_electrodes]);
reps_to_use = 1:n_reps_total;
% reps_to_use = reps_to_use(randperm(n_reps_total));
reps_to_use = reps_to_use(1:I.n_reps);
raster_rep1 = squeeze_dims(nanmean(raster_separate_reps(:,reps_to_use(1:2:end),:,:),2),2);
raster_rep2 = squeeze_dims(nanmean(raster_separate_reps(:,reps_to_use(2:2:end),:,:),2),2);

% test retest correlation
r = fastcorr(reshape(raster_rep1, [n_bins * n_conditions, n_electrodes]), ...
    reshape(raster_rep2, [n_bins * n_conditions, n_electrodes]));

% plot
if I.plot
    figure;
    bar(r, 'FaceColor', [1 1 1]*0.5)
    xlim([0, n_electrodes+1]);
    xlabel('Units');
    ylabel('Correlation (r)');
    title(sprintf('reliability based on %d repetitions\n', I.n_reps));
    
    % save
    if ~isempty(I.figure_directory)
        fname = [I.figure_directory '/' recording_id '-thresh'  ...
            num2str(I.spike_thresh) '-' num2str(I.fwhm_ms) 'ms'  ...
            '-raster' num2str(I.raster_sr) 'Hz-reps' num2str(I.n_reps)];
        if I.single_unit
            fname = [fname '-sorted'];
        end
        if exist('export_fig', 'file')
            export_fig([fname '.pdf'], '-pdf', '-transparent');
            export_fig([fname '.png'], '-png', '-transparent', '-r100');
        end
    end
end

%% Plot individual electrodes

if ~isempty(I.individ_units)
    figh = figure;
    for j = 1:length(I.individ_units)
        
        unit = I.individ_units(j);
        Y1 = raster_rep1(:,:,unit);
        Y2 = raster_rep2(:,:,unit);
        
        % plot
        clf(figh);
        subplot(2,1,1);
        imagesc(Y1', [-1 1] * quantile(raster_smooth(:), 0.99))
        subplot(2,1,2);
        imagesc(Y2', [-1 1] * quantile(raster_smooth(:), 0.99))
        for i = 1:2
            subplot(2,1,i);
            colormap(flipud(cbrewer('div', 'RdBu', 128)));
            xticks = get(gca, 'XTick');
            set(gca, 'XTick', xticks, 'XTickLabel', tbins(xticks));
            ylabel('Stimuli'); xlabel('Time (s)');
            title(sprintf('unit %d, split %d', unit, i));
        end
        
        % save
        if ~isempty(I.figure_directory)
            test_retest_individ_elec_directory = [I.figure_directory ...
                '/rasters-' recording_id '-I.thresh'  ...
                num2str(I.spike_thresh) '-' num2str(I.fwhm_ms) 'ms' ...
                '-raster' num2str(I.raster_sr) 'Hz-reps' num2str(I.n_reps)];
            if ~exist(test_retest_individ_elec_directory, 'dir'); mkdir(test_retest_individ_elec_directory); end
            if I.single_unit
                fname = [test_retest_individ_elec_directory '/su' num2str(unit)];
            else
                fname = [test_retest_individ_elec_directory '/elec' num2str(unit)];
            end
            export_fig([fname '.pdf'], '-pdf', '-transparent');
            export_fig([fname '.png'], '-png', '-transparent', '-r100');
        end
    end
end
