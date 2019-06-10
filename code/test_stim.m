% Some code to bug-check the stimuli

elec = 1;

prestim_silence = 0.5;

t_old = (0:length(X.quilt)-1)'/40e3;
t_new = (0:1/raster_sr:10-1/raster_sr)';
raster_smooth = nan(length(t_new), n_trials);


for i = 1:n_trials
    X = load(['/Users/svnh2/Dropbox (MIT)/mindhive/scrambling-ferrets/stimuli' ...
        '/naturalsound20-quilt-0.5sec-catmethod2' ...
        '/seg-' num2str(round(segs(i))) 'ms-stim' ...
        num2str(stims(i)) '-order' num2str(orders(i)) '.mat']);
    
    y = X.quilt;
    [B,A] = butter(4, [4000, 4300]/(20e3)); 
    y = filtfilt(B, A, y);
    y = abs(hilbert(y)).^0.3;
    raster_smooth(:,i) = resample(y, raster_sr, 40e3);%interp1(t_old, X.quilt, t_new);
end

%%
Z = raster_smooth(:,trial_order(1:8:end));
plot(Z)

nfft = 2^9;
F = nan(nfft/2+1, size(Z,2));
for i = 1:12
    [F(:,i), f] = fftplot2(Z(:,i), raster_sr, 'nfft', nfft);
end

F = F(2:end, :);
f = f(2:end);

Fdb = 10*log10(bsxfun(@times, F, f'));

imagesc(Fdb);

% imagesc(raster_smooth)

%%
raster_smooth = [zeros(prestim_silence*raster_sr, n_trials); raster_smooth; zeros(prestim_silence*raster_sr, n_trials)];
bins = [ flipud(-(1:prestim_silence*raster_sr)'/raster_sr); t_new; t_new(end) + (1:prestim_silence*raster_sr)'/raster_sr];


%%


%%


fname1 = ['/Users/svnh2/Dropbox (MIT)/mindhive/scrambling-ferrets/stimuli' ...
        '/naturalsound20-quilt-0.5sec-catmethod2' ...
        '/seg-500ms-stim1-order1.mat'];
X1 = load(fname1);

fname2 = ['/Users/svnh2/Dropbox (MIT)/mindhive/scrambling-ferrets/stimuli' ...
        '/naturalsound20-quilt-0.5sec-catmethod2' ...
        '/seg-500ms-stim1-order2.mat'];
X2 = load(fname2);


clc;
[X1.order, X2.order]

plot(bins, [y1,y2]);
y1 = mean(raster_smooth(:,segs==500 & orders==1),2);
y2 = mean(raster_smooth(:,segs==500 & orders==2),2);



% sound(X2.quilt(1:40e3*2), 40e3)
    
    


%%


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
            Y1(k, :, :) = interp1( bins', X1, lag_t' + seg_dur_sec*(k-1));
            Y2(k, :, :) = interp1( bins', X2, lag_t' + seg_dur_sec*(k-1));
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
    '-thresh' num2str(thresh) '-' num2str(fwhm_ms) 'ms' '-raster' num2str(raster_sr) 'Hz'];
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



%%




%%
r_same_order = nan(n_lags, n_stims, n_seg);
for i = 1:n_seg
    for j = 1:n_stims
        
        % responses for first order averaged separately across even and odd reps
        xi = segs == round(unique_segs(i));
        
        
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
