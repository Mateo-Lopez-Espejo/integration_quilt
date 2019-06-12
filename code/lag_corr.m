function lag_corr(recording_directory, recording_id, varargin)

% Primary analysis script for performing the lagged correlation analysis.

% 2018-03-26: Created, Sam NH
 
%% Parameters

I.figure_directory = '';

I.lag_win = [0, 1];

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

if isempty(I.units)
    I.units = 1:size(raster_smooth,3);
end

%% Do lag analysis

figh = figure;
for elec = I.units(:)'
    
    %% Measured lagged correlation
    
    % total duration of the stimulus in seconds
    stim_dur = 10;
    
    % segment durations
    unique_segs = [15.625, 31.25, 62.50, 125, 250, 500];
    assert(all(round(unique_segs) == T.unique_segs));
    
    % window over which to compute lags
    lag_t = I.lag_win(1):1/I.raster_sr:I.lag_win(2);
    n_lags = length(lag_t);
    
    r_same_order = nan(n_lags, T.n_stims, T.n_seg);
    r_diff_order = nan(n_lags, T.n_stims, T.n_seg);
    r_diff_order_matched = nan(n_lags, T.n_stims, T.n_seg);
    for i = 1:T.n_seg
        for j = 1:T.n_stims
            
            % responses for first order averaged separately across even and odd reps
            xi = T.segs == round(unique_segs(i)) & T.stims == T.unique_stims(j) & T.orders == 1;
            X1 = raster_smooth(:, xi, elec);
            X1 = cat(2, nanmean(X1(:,1:2:end,:),2), nanmean(X1(:,2:2:end,:),2));
            
            % responses for second order averaged separately across even and odd reps
            xi = T.segs == round(unique_segs(i)) & T.stims == T.unique_stims(j) & T.orders == 2;
            X2 = raster_smooth(:, xi, elec);
            X2 = cat(2, nanmean(X2(:,1:2:end,:),2), nanmean(X2(:,2:2:end,:),2));
            
            % load the segment T.orders
            Z1 = load(['/Users/svnh2/Dropbox (MIT)/mindhive/scrambling-ferrets/stimuli' ...
                '/naturalsound20-quilt-0.5sec-catmethod2' ...
                '/seg-' num2str(round(unique_segs(i))) 'ms-stim' num2str(T.unique_stims(j)) '-order1.mat']);
            Z2 = load(['/Users/svnh2/Dropbox (MIT)/mindhive/scrambling-ferrets/stimuli' ...
                '/naturalsound20-quilt-0.5sec-catmethod2' ...
                '/seg-' num2str(round(unique_segs(i))) 'ms-stim' num2str(T.unique_stims(j)) '-order2.mat']);
            
            % number of segments per stim and duration of the segment in seconds
            T.n_seg_per_stim = stim_dur/(unique_segs(i)/1000);
            seg_dur_sec = (unique_segs(i)/1000);
            
            % compute the response to each segment
            Y1 = nan(T.n_seg_per_stim, n_lags, 2);
            Y2 = nan(T.n_seg_per_stim, n_lags, 2);
            for k = 1:T.n_seg_per_stim
                Y1(k, :, :) = interp1( tbins' - T.prestim_silence, X1, lag_t' + seg_dur_sec*(k-1));
                Y2(k, :, :) = interp1( tbins' - T.prestim_silence, X2, lag_t' + seg_dur_sec*(k-1));
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
        
    % plot correlaiton vs. lag
    clf(figh);
    set(figh, 'Position', [100, 100, 500, 900]);
    corr_range = [0 1];
    for k = 1:length(unique_segs)
        subplot(length(unique_segs), 1, k);
        plot(lag_t * 1000, [r_same_order(:,k), r_diff_order(:,k), r_diff_order_matched(:,k)])
        hold on;
        plot(unique_segs(k)*[1 1], corr_range, 'k--');
        xlim(I.lag_win * 1000);
        ylim(corr_range);
        xlabel('Time Lag (ms)');
        ylabel('Pearson Correlation');
        title(sprintf('%.0f ms', unique_segs(k)))
    end
    
    if ~isempty(I.figure_directory)
        lag_directory = [I.figure_directory '/' recording_id ...
            '-I.thresh' num2str(I.spike_thresh) '-' num2str(I.fwhm_ms) 'ms' '-raster' num2str(I.raster_sr) 'Hz'];
        if ~exist(lag_directory, 'dir'); mkdir(lag_directory); end
        if I.single_unit
            export_fig([lag_directory '/su' num2str(elec) '.pdf'], '-pdf', '-transparent');
        else
            export_fig([lag_directory '/elec' num2str(elec) '.pdf'], '-pdf', '-transparent');
        end
    end
    
    %% Normalized similarity measure
    
    % measure to predict with the model
    z_same_order = atanh(r_same_order);
    z_diff_order_matched = atanh(r_diff_order_matched);
    Y = nan(size(z_same_order));
    for k = 1:length(unique_segs)
        baseline = mean(z_same_order(:,k));
        Y(:,k) = 1 - (z_same_order(:,k)-z_diff_order_matched(:,k))/baseline;
    end
    
    %% Model prediction
    
    % segment dependent weights
    w = sqrt([1280, 640, 320, 160, 80, 40]-3);
    w = w / sum(w);
    
    % rf_dur_sec = unique_segs/1000;
    rf_dur_sec = logspace(log10(unique_segs(1)/1000), log10(unique_segs(end)/1000), 20);
    delay_dur_smp = 1:0.5*I.raster_sr;
    err = nan(length(rf_dur_sec), length(delay_dur_smp));
    Yh = nan([size(Y),length(rf_dur_sec), length(delay_dur_smp)]);
    for i = 1:length(rf_dur_sec)
        
        fprintf('%d\n', i); drawnow;
        for j = 1:length(delay_dur_smp)
            for k = 1:length(unique_segs)
                
                % time / sample vector
                t = -2*I.lag_win(2)*I.raster_sr:2*I.lag_win(2)*I.raster_sr;
                
                % gaussian window
                % convert to sigma (equivalent rectangular window)
                sig_sec = rf_dur_sec(i) / 2.5066;
                sig_smp = I.raster_sr * sig_sec;
                % sig = 1;
                % equiv_rec_width = 1/normpdf(0, 0, sig)
                h = normpdf(t, 0, sig_smp);
                h = h / sum(h);
                % plot(t, h)
                
                % rectangular segment
                seg = zeros(size(t));
                seg_dur_smp = round(I.raster_sr * unique_segs(k)/1000);
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
    
    %% Plot model stats
    
    clf(figh);
    set(figh, 'Position', [100, 100, 900, 900]);
    corr_range = [-0.5 1.5];
    for k = 1:length(unique_segs)
        subplot(4, 2, k);
        plot(I.lag_win * 1000, [1,1], 'k--', 'LineWidth', 2);
        hold on;
        plot(I.lag_win * 1000, [0,0], 'k--', 'LineWidth', 2);
        h1 = plot(lag_t * 1000, Y(:,k), 'LineWidth', 2);
        h2 = plot(lag_t * 1000, Ybest(:,k), 'LineWidth', 2);
        plot(unique_segs(k)*[1 1], corr_range, 'k--', 'LineWidth', 2);
        xlim(I.lag_win * 1000);
        ylim(corr_range);
        xlabel('Time Lag (ms)');
        ylabel('Pearson Correlation');
        title(sprintf('Seg: %.0f ms', unique_segs(k)))
        if k == 1
            legend([h1, h2], 'Data', 'Model');
        end
        box off;
    end
    if ~isempty(I.figure_directory)
        if I.single_unit
            fname = [lag_directory '/su' num2str(elec) '-model-prediction-lineplot'];
        else
            fname = [lag_directory '/elec' num2str(elec) '-model-prediction-lineplot'];
        end
        export_fig([fname '.pdf'], '-pdf', '-transparent');
    end
    
    % plot prediction for best delay
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
            title(sprintf('rf=%.f ms, delay=%.f ms', best_rf_sec*1000, best_delay_smp/I.raster_sr*1000));
        end
    end
    if ~isempty(I.figure_directory)
        if I.single_unit
            fname = [lag_directory '/su' num2str(elec) '-model-prediction'];
        else
            fname = [lag_directory '/elec' num2str(elec) '-model-prediction'];
        end
        export_fig([fname '.png'], '-png', '-transparent', '-r100');
    end

    
    % plot the error
    clf(figh);
    set(figh, 'Position', [100, 100, 600, 600]);
    imagesc(-err, quantile(-err(:), [0.8, 1]));
    colormap(parula(128));
    colorbar;
    ytick = interp1(log2(rf_dur_sec), 1:length(rf_dur_sec), log2(rf_dur_sec));
    set(gca, 'YTick', ytick, 'YTickLabel', rf_dur_sec*1000);
    xticks = get(gca, 'XTick');
    set(gca, 'XTickLabel', 1000*delay_dur_smp(xticks)/I.raster_sr);
    xlabel('Delay (ms)'); ylabel('Receptive Field Duration (ms)');
    title(sprintf('rf=%.f ms, delay=%.f ms', best_rf_sec*1000, best_delay_smp/I.raster_sr*1000));
    set(gca, 'FontSize', 12);
    if ~isempty(I.figure_directory)
        if I.single_unit
            fname = [lag_directory '/su' num2str(elec) '-model-error'];
        else
            fname = [lag_directory '/elec' num2str(elec) '-model-error'];
        end
        export_fig([fname '.png'], '-png', '-transparent', '-r100');
    end

    
end
