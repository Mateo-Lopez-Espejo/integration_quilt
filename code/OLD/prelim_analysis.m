% function prelim_analysis(path, 

recording_id = 'tomette014a09_p_NSD';
recording_directory = '/Volumes/data/Tomette/tomette014';

% Preliminary analysis of the scrambling experiment data
% 
% 2017-12-05: Created Sam NH

% paths to external code repositories
addpath(genpath([root_directory '/general-analysis-code']));
addpath(genpath([root_directory '/export_fig_v3']));

project_directory = [root_directory '/scrambling-ferrets'];

%% Plot PSTH to detect onset responses

figure_directory = [project_directory '/figures/PSTH'];
if ~exist(figure_directory, 'dir'); mkdir(figure_directory); end

% recordings
% recording_ids = {...
%     'pelardon060a03_p_NSD', 'pelardon062a01_p_NSD', 'pelardon063a02_p_NSD', ...
%     'pelardon064a01_p_NSD', 'pelardon065a05_p_NSD', 'pelardon066a02_p_NSD'};
% recording_ids = {'boulette030a10_p_NSD'};
% recording_ids = {'tomette002a10_p_NSD'};
% recording_ids = {'tomette010a08_p_NSD'};
% recording_ids = {'tomette012a06_p_NSD'};
% recording_ids = {'tomette014a09_p_NSD'};
% animal = 'Tomette';
% folder = 'tomette014';
    
% sampling rate of the spikes
% spike_sr = 31250;

% number of electrodes in the array
% n_electrodes = 32;

% threshold
I.thresh = -4;

%% Parse stimulus info

run([recording_directory '/' recording_id '.m']);

% total numer of trials
n_trials = exptevents(end).Trial;
% n_trials = 132;

% exptevents struct into stimulus labels
notes = {exptevents(:).Note};
start_times = cat(1, exptevents(:).StartTime);
trial_indices = cat(1, exptevents(:).Trial);
stim_labels = cell(1, n_trials);
onset_times = nan(1, n_trials);
for i = 1:n_trials
    xi = find(trial_indices == i);
    for j = 1:length(xi)
        s = regexp(notes{xi(j)}, 'seg[-\w]*', 'match');
        if ~isempty(s)
            stim_labels(i) = s;
        end
        if ~isempty(s) && strcmp(notes{xi(j)}(1:4), 'Stim')
            onset_times(i) = start_times(xi(j));
        end
    end
end

% parse stimulus labels into durations / stims / orders
segs = nan(1, n_trials);
stims = nan(1, n_trials);
orders = nan(1, n_trials);
for i = 1:n_trials
    x = regexp(stim_labels{i}, '-', 'split');
    segs(i) = round(str2double(x{2}(1:end-2)));
    stims(i) = round(str2double(x{3}(end)));
    orders(i) = round(str2double(x{4}(end)));
end

% come up with a trial order based on the above info
trial_order = nan(1, n_trials);
unique_segs = sort(unique(segs));
n_seg = length(unique_segs);
unique_stims = sort(unique(stims));
n_stims = length(unique(stims));
unique_orders = sort(unique(orders));
n_orders = length(unique(orders));
count = 0;
for i = 1:n_seg
    xi = find(unique_segs(i) == segs);
    for j = 1:n_stims
        yi = find(stims(xi) == unique_stims(j));
        for k = 1:n_orders
            zi = find(orders(xi(yi)) == unique_orders(k));
            trial_order((1:length(zi)) + count) = xi(yi(zi));
            count = count + length(zi);
        end
    end    
end

prestim_silence = onset_times(1);
assert(all(abs(onset_times - prestim_silence) < 1e-6))

%% PSTH

% window relative to trial onset (not stim onset)
psth_win = [0, 3];
psth_sr = 50;

figh = figure;
set(figh, 'Position', [100 100 1400 800]);
for i = 1%:length(recording_ids)
    
    recording_id = recording_ids{i};
    
    clf(figh);
    for elec = 1:n_electrodes
        
        % load data for this electrode
        X = load(['/Volumes/data/' animal '/' folder '/tmp/' ...
            recording_id '.001.1.elec' num2str(elec) '.sig' num2str(I.thresh) '.NOCOM.mat']);
        
        % calculate psth
        bins = psth_win(1):1/psth_sr:psth_win(2);
        [N, p] = myhist(X.spikebin/spike_sr, bins);
        
        % plot
        n_cols = round(sqrt(n_electrodes));
        n_rows = ceil(n_electrodes/n_cols);
        subplot(n_rows, n_cols, elec);
        plot(bins, p, 'k-', 'LineWidth', 2);
        xlim(psth_win);
        ylim([0, max(p)]);
        title(sprintf('electrode %d', elec));
        xlabel('Time (sec)');
        ylabel('Probability of spike');
        hold on;
        plot(prestim_silence*[1, 1], [0, max(p)], 'r--', 'LineWidth', 2);
        box off;
        
        
    end
    
    % save
    export_fig([figure_directory '/' recording_ids{i} '-I.thresh'  num2str(I.thresh) '.pdf'], '-pdf', '-transparent');
    
end

%% Time analysis

% window to analyze
win = [0, 10+prestim_silence*2+0.5];

% sampling rate of the raster
raster_sr = 500;

% bins for the raster
bins = win(1) : 1/raster_sr : win(2);
n_bins = length(bins);

% raster for all units
raster = nan(n_bins, n_trials, n_electrodes);
for elec = 1:n_electrodes
    
    % load data
    load(['/Volumes/data/' animal '/' folder '/tmp' '/' ...
        recording_id '.001.1.elec' num2str(elec) '.sig' num2str(I.thresh) '.NOCOM.mat'], ...
        'trialid', 'spikebin');
    
    % compute raster
    for i = 1:n_trials
        spike_times = spikebin(trialid == i)/spike_sr;
        if ~isempty(spike_times)
            raster(:, i, elec) = myhist(spike_times, bins);
        end
    end
end

% smooth with Gaussian kernel
fwhm_ms = 10;
if fwhm_ms > 0
    raster_smooth = mysmooth(raster, raster_sr, fwhm_ms);
else
    raster_smooth = raster;
end

%% Plot test retest reliability

n_conditions = n_seg * n_orders * n_stims;
n_reps = n_trials/n_conditions;

% separate out different trials
% time/bins x repetitions x conditions x electrodes
raster_shape = reshape(raster_smooth(:,trial_order,:), [n_bins, n_reps, n_conditions, n_electrodes]);

% separ
raster1 = squeeze_dims(nanmean(raster_shape(:,1:2:end,:,:),2),2);
raster2 = squeeze_dims(nanmean(raster_shape(:,2:2:end,:,:),2),2);

% test retets correlation
r = fastcorr(reshape(raster1, [n_bins * n_conditions, n_electrodes]), ...
    reshape(raster2, [n_bins * n_conditions, n_electrodes]));

% plot
figure;
bar(r, 'FaceColor', [1 1 1]*0.5)
xlim([0, n_electrodes+1]);
xlabel('Electrodes');
ylabel('Correlation (r)');

% save
test_retest_directory = [project_directory '/figures/test-retest'];
if ~exist(test_retest_directory, 'dir'); mkdir(test_retest_directory); end
fname = [test_retest_directory '/' recording_id '-I.thresh'  ...
    num2str(I.thresh) '-' num2str(fwhm_ms) 'ms' '-raster' num2str(raster_sr) 'Hz'];
export_fig([fname '.pdf'], '-pdf', '-transparent');
export_fig([fname '.png'], '-png', '-transparent', '-r100');

%% Plot the reliability for individual electrodes

figh = figure;
for elec = [18 12 10 25]%[23,10,14,16]
    Y1 = raster1(:,:,elec);
    Y2 = raster2(:,:,elec);
    
    % plot
    clf(figh);
    bounds = [-1 1] * quantile([Y1(:); Y2(:)], 0.99);
    subplot(2,1,1);
    imagesc(Y1', [-1 1] * quantile(raster_smooth(:), 0.99))
    subplot(2,1,2);
    imagesc(Y2', [-1 1] * quantile(raster_smooth(:), 0.99))
    for i = 1:2
        subplot(2,1,i);
        colormap(flipud(cbrewer('div', 'RdBu', 128)));
        xticks = get(gca, 'XTick');
        set(gca, 'XTick', xticks, 'XTickLabel', bins(xticks));
        ylabel('Stimuli'); xlabel('Time (s)');
        title(sprintf('elec %d, split %d', elec, i));
    end
    
    % save
    test_retest_individ_elec_directory = [project_directory '/figures' ...
        '/raster-test-retest/' recording_id '-I.thresh'  ...
        num2str(I.thresh) '-' num2str(fwhm_ms) 'ms' '-raster' num2str(raster_sr) 'Hz'];
    if ~exist(test_retest_individ_elec_directory, 'dir'); mkdir(test_retest_individ_elec_directory); end
    fname = [test_retest_individ_elec_directory '/elec' num2str(elec)];
    export_fig([fname '.pdf'], '-pdf', '-transparent');
    export_fig([fname '.png'], '-png', '-transparent', '-r100');
end

%% segment

for elec = [18 12 10 25]
    
    % total duration of the stimulus in seconds
    stim_dur = 10;
    
    % segment durations
    unique_segs = [15.625, 31.25, 62.50, 125, 250, 500];
    
    % window over which to compute lags
    lag_win = [0, 1];
    lag_t = lag_win(1):1/raster_sr:lag_win(2);
    n_lags = length(lag_t);
    
    r_same_order = nan(n_lags, n_stims, n_seg);
    r_diff_order = nan(n_lags, n_stims, n_seg);
    r_diff_order_matched = nan(n_lags, n_stims, n_seg);
    for i = 1:n_seg
        for j = 1:n_stims
            
            % responses for first order averaged separately across even and odd reps
            xi = segs == round(unique_segs(i)) & stims == unique_stims(j) & orders == 1;
            X1 = raster_smooth(:, xi, elec);
            X1 = cat(2, nanmean(X1(:,1:2:end,:),2), nanmean(X1(:,2:2:end,:),2));
            
            % responses for second order averaged separately across even and odd reps
            xi = segs == round(unique_segs(i)) & stims == unique_stims(j) & orders == 2;
            X2 = raster_smooth(:, xi, elec);
            X2 = cat(2, nanmean(X2(:,1:2:end,:),2), nanmean(X2(:,2:2:end,:),2));
            
            % load the segment orders
            Z1 = load(['/Users/svnh2/Dropbox (MIT)/mindhive/scrambling-ferrets/stimuli' ...
                '/naturalsound20-quilt-0.5sec-catmethod2' ...
                '/seg-' num2str(round(unique_segs(i))) 'ms-stim' num2str(unique_stims(j)) '-order1.mat']);
            Z2 = load(['/Users/svnh2/Dropbox (MIT)/mindhive/scrambling-ferrets/stimuli' ...
                '/naturalsound20-quilt-0.5sec-catmethod2' ...
                '/seg-' num2str(round(unique_segs(i))) 'ms-stim' num2str(unique_stims(j)) '-order2.mat']);
            
            % number of segments per stim and duration of the segment in seconds
            n_seg_per_stim = stim_dur/(unique_segs(i)/1000);
            seg_dur_sec = (unique_segs(i)/1000);
            
            % compute the response to each segment
            Y1 = nan(n_seg_per_stim, n_lags, 2);
            Y2 = nan(n_seg_per_stim, n_lags, 2);
            for k = 1:n_seg_per_stim
                Y1(k, :, :) = interp1( bins' - prestim_silence, X1, lag_t' + seg_dur_sec*(k-1));
                Y2(k, :, :) = interp1( bins' - prestim_silence, X2, lag_t' + seg_dur_sec*(k-1));
            end
            
            % correlation for same vs. different order
            r_same_order(:,j,i) = nanfastcorr(Y1(:,:,1), Y1(:,:,2))/2 + nanfastcorr(Y2(:,:,1), Y2(:,:,2))/2;
            r_diff_order(:,j,i) ...
                = nanfastcorr(Y1(:,:,1), Y2(:,:,1))/4 ...
                + nanfastcorr(Y1(:,:,1), Y2(:,:,2))/4 ...
                + nanfastcorr(Y1(:,:,2), Y2(:,:,1))/4 ...
                + nanfastcorr(Y1(:,:,2), Y2(:,:,2))/4;
            
            % sort by the segment order
            [~,xi1] = sort(Z1.order);
            [~,xi2] = sort(Z2.order);
            assert(all(Z1.order(xi1) == Z2.order(xi2)));
            
            % recompute after matching segments
            r_diff_order_matched(:,j,i) ...
                = nanfastcorr(Y1(xi1,:,1), Y2(xi2,:,1))/4 ...
                + nanfastcorr(Y1(xi1,:,1), Y2(xi2,:,2))/4 ...
                + nanfastcorr(Y1(xi1,:,2), Y2(xi2,:,1))/4 ...
                + nanfastcorr(Y1(xi1,:,2), Y2(xi2,:,2))/4;
            
        end
    end
    
    % average across the two stimuli
    r_same_order = squeeze_dims(nanmean(r_same_order,2),2);
    r_diff_order = squeeze_dims(nanmean(r_diff_order,2),2);
    r_diff_order_matched = squeeze_dims(nanmean(r_diff_order_matched,2),2);
    
    lag_directory = [project_directory '/figures/lag-correlation/' recording_id ...
        '-I.thresh' num2str(I.thresh) '-' num2str(fwhm_ms) 'ms' '-raster' num2str(raster_sr) 'Hz'];
    if ~exist(lag_directory, 'dir'); mkdir(lag_directory); end
    
    % plot correlaiton vs. lag
    figh = figure;
    set(figh, 'Position', [100, 100, 500, 900]);
    corr_range = [0 1];
    for k = 1:length(unique_segs)
        subplot(length(unique_segs), 1, k);
        plot(lag_t * 1000, [r_same_order(:,k), r_diff_order(:,k), r_diff_order_matched(:,k)])
        hold on;
        plot(unique_segs(k)*[1 1], corr_range, 'k--');
        xlim(lag_win * 1000);
        ylim(corr_range);
        xlabel('Time Lag (ms)');
        ylabel('Pearson Correlation');
        title(sprintf('%.0f ms', unique_segs(k)))
        
    end
    
    export_fig([lag_directory '/elec' num2str(elec) '.pdf'], '-pdf', '-transparent')
    export_fig([lag_directory '/elec' num2str(elec) '.png'], '-png', '-transparent', '-r100')
    
    %%
    
    % measure to predict with the model
    z_same_order = atanh(r_same_order);
    z_diff_order_matched = atanh(r_diff_order_matched);
    Y = nan(size(z_same_order));
    for k = 1:length(unique_segs)
        baseline = mean(z_same_order(:,k));
        Y(:,k) = 1 - (z_same_order(:,k)-z_diff_order_matched(:,k))/baseline;
    end
    
    
    %%
    
    % segment dependent weights
    w = sqrt([1280, 640, 320, 160, 80, 40]-3);
    w = w / sum(w);
    
    % rf_dur_sec = unique_segs/1000;
    rf_dur_sec = logspace(log10(unique_segs(1)/1000), log10(unique_segs(end)/1000), 20);
    delay_dur_smp = 1:0.5*raster_sr;
    err = nan(length(rf_dur_sec), length(delay_dur_smp));
    Yh = nan([size(Y),length(rf_dur_sec), length(delay_dur_smp)]);
    for i = 1:length(rf_dur_sec)
        
        fprintf('%d\n', i); drawnow;
        for j = 1:length(delay_dur_smp)
            for k = 1:length(unique_segs)
                
                % time / sample vector
                t = -2*lag_win(2)*raster_sr:2*lag_win(2)*raster_sr;
                
                % gaussian window
                % convert to sigma (equivalent rectangular window)
                sig_sec = rf_dur_sec(i) / 2.5066;
                sig_smp = raster_sr * sig_sec;
                % sig = 1;
                % equiv_rec_width = 1/normpdf(0, 0, sig)
                h = normpdf(t, 0, sig_smp);
                h = h / sum(h);
                % plot(t, h)
                
                % rectangular segment
                seg = zeros(size(t));
                seg_dur_smp = round(raster_sr * unique_segs(k)/1000);
                seg(t <= seg_dur_smp & t >= 0) = 1;
                % plot(t, [h'/max(h(:)), seg'/max(seg(:))]);
                
                % convolve receptive field with segment
                area = myconv(seg', h', 'causal', false);
                % plot(t, area);
                % ylim([0 1]);
                
                predictor = area(t >= -delay_dur_smp(j));
                Yh(:,k,i,j) = predictor(1:size(Yh,1));
                
            end
            E = Yh(:,:,i,j) - Y;
            Ew = bsxfun(@times, E, w);
            err(i,j) = mean(Ew(:).^2);
        end
    end
    
    % find best prediction
    [~,xi] = min(err(:));
    [a,b] = ind2sub(size(err), xi);
    Ybest = Yh(:,:,a,b);
    best_rf_sec = rf_dur_sec(a);
    best_delay_smp = delay_dur_smp(b);
    
    %%
    
    figh = figure;
    clf(figh);
    set(figh, 'Position', [100, 100, 900, 900]);
    corr_range = [-0.5 1.5];
    for k = 1:length(unique_segs)
        subplot(4, 2, k);
        plot(lag_win * 1000, [1,1], 'k--', 'LineWidth', 2);
        hold on;
        plot(lag_win * 1000, [0,0], 'k--', 'LineWidth', 2);
        h1 = plot(lag_t * 1000, Y(:,k), 'LineWidth', 2);
        h2 = plot(lag_t * 1000, Ybest(:,k), 'LineWidth', 2);
        plot(unique_segs(k)*[1 1], corr_range, 'k--', 'LineWidth', 2);
        xlim(lag_win * 1000);
        ylim(corr_range);
        xlabel('Time Lag (ms)');
        ylabel('Pearson Correlation');
        title(sprintf('Seg: %.0f ms', unique_segs(k)))
        if k == 1
            legend([h1, h2], 'Data', 'Model');
        end
        box off;
    end
    fname = [lag_directory '/elec' num2str(elec) '-model-prediction-lineplot'];
    export_fig([fname '.png'], '-png', '-transparent', '-r100');
    export_fig([fname '.pdf'], '-pdf', '-transparent');
    
    % plot prediction for best delay
    figh = figure;
    clf(figh);
    set(figh, 'Position', [100, 100, 900, 600]);
    subplot(2,1,1);
    imagesc(Y', [-1 1]);
    subplot(2,1,2);
    imagesc(Ybest', [-1 1]);
    for i = 1:2
        subplot(2,1,i);
        colorbar;
        colormap(flipud(cbrewer('div', 'RdBu', 128)));
        set(gca, 'YTick', 1:length(unique_segs), 'YTickLabel', unique_segs);
        xticks = get(gca, 'XTick');
        set(gca, 'XTick', xticks+1, 'XTickLabel', lag_t(xticks+1)*1000);
        ylabel('Seg Dur (ms)'); xlabel('Lag (ms)');
        if i == 2
            title(sprintf('rf=%.f ms, delay=%.f ms', best_rf_sec*1000, best_delay_smp/raster_sr*1000));
        end
    end
    fname = [lag_directory '/elec' num2str(elec) '-model-prediction'];
    export_fig([fname '.png'], '-png', '-transparent', '-r100');
    
    % plot the error
    figh = figure;
    clf(figh);
    set(figh, 'Position', [100, 100, 600, 600]);
    imagesc(-err, quantile(-err(:), [0.8, 1]));
    colormap(parula(128));
    colorbar;
    ytick = interp1(log2(rf_dur_sec), 1:length(rf_dur_sec), log2(rf_dur_sec));
    set(gca, 'YTick', ytick, 'YTickLabel', rf_dur_sec*1000);
    xticks = get(gca, 'XTick');
    set(gca, 'XTickLabel', 1000*delay_dur_smp(xticks)/raster_sr);
    xlabel('Delay (ms)'); ylabel('Receptive Field Duration (ms)');
    title(sprintf('rf=%.f ms, delay=%.f ms', best_rf_sec*1000, best_delay_smp/raster_sr*1000));
    set(gca, 'FontSize', 12);
    fname = [lag_directory '/elec' num2str(elec) '-model-error'];
    export_fig([fname '.png'], '-png', '-transparent', '-r100');
    
end


%% Onset responses


% total duration of the stimulus in seconds
stim_dur = 10;

% segment durations
unique_segs = [15.625, 31.25, 62.50, 125, 250, 500];

% window over which to compute lags
lag_win = [-0.5, 0.5];
lag_t = lag_win(1):1/raster_sr:lag_win(2);
n_lags = length(lag_t);

y = nan(n_lags-1, n_seg);

for i = 1:n_seg
    
    
    % responses for first order averaged separately across even and odd reps
    
    elec = 23;
    
    xi = segs == round(unique_segs(i));
    X1 = raster_smooth(:, xi, elec);
    
    % number of segments per stim and duration of the segment in seconds
    n_seg_per_stim = stim_dur/(unique_segs(i)/1000);
    seg_dur_sec = (unique_segs(i)/1000);
    
    
    Y = nan(n_lags, n_seg_per_stim);
    for k = 1:n_seg_per_stim
        Y(:,k) = mean(interp1( bins' - prestim_silence, X1, lag_t' + seg_dur_sec*(k-1)),2);
    end
    
    y(:,i) = nanmean(abs(diff(Y,[],1)),2);
end


figure;
plot(lag_t(1:n_lags-1), y)


