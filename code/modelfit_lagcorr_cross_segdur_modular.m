function M = modelfit_lagcorr_cross_segdur_modular(L, varargin)

% Fit model to lagged correlation data to quantitatively estimate 
% integration window
% 
% 2019-04-05: Last edited, Sam NH

% optional parameters
I.plot_win = L.lag_t([1 end]);
I.plot_figure = true;
I.keyboard = false;
I.overwrite = false;
I.normcorr = false;
I.divnorm = false;
I.mindenom = 0.01;
I.weightdenom = false;
I.distr = 'gamma';
I.tranweightnsegs = 'none'; % applied to total segs
I.tranweightdenom = 'none'; % applied to denominator
I.forcecausal = false;
I.bestcausal = true; % select only amongst causal solutions (not relevant for Gausssian)
I.plotcausal = true; % only plot errors for causal solutions
switch I.distr
    case 'gauss'
        I.shape = 1;
    case 'gamma'
        I.shape = logspace(log10(1), log10(10), 10);
    otherwise
        error('No matching distribution');
end
[I, ~, C_value] = parse_optInputs_keyvalue(varargin, I);

if I.keyboard
    keyboard;
end

% string identifying parameters of model fit
param_string_modelfit = [...
    struct2string(I, 'include_fields', {'distr', 'normcorr'}), ...
    '_' struct2string(C_value, 'include_fields', ...
    {'divnorm', 'mindenom', 'weightdenom','tranweightnsegs','tranweightdenom','shape'})];
if param_string_modelfit(end)=='_'; param_string_modelfit = param_string_modelfit(1:end-1); end

% function used to transform correlation
switch I.tranweightnsegs
    case 'none'
        tranweightnsegsfn = @(x)x;
    case 'sqrt'
        tranweightnsegsfn = @(x)sqrt(x);
    case 'sqrtm3'
        tranweightnsegsfn = @(x)sqrt(x-3);
    otherwise
        error('No matching tranweightnsegs parameter');
end

% function used to transform weights
switch I.tranweightdenom
    case 'none'
        tranweightdenomfn = @(x)x;
    case 'square'
        tranweightdenomfn = @(x)x.^2;
    otherwise
        error('No matching tranweightdenom parameter');
end

%% Model fit

MAT_file = [L.output_directory '/model_fit_' param_string_modelfit '.mat'];
if ~exist(MAT_file, 'file') || I.overwrite
    
    [n_lags, n_seg, n_channels] = size(L.same_context);
    assert(n_channels == length(L.channels));
    assert(length(L.unique_segs)==n_seg);
    assert(n_lags == length(L.lag_t))
    
    % measure to predict with the model
    same_context = L.same_context;
    diff_context = L.diff_context;
    denom_factor = nan(size(L.same_context));
    M.Y = nan(size(same_context));
    for k = 1:n_seg
        if I.normcorr
            assert(~I.divnorm);
            baseline = mean(same_context(:,k,:));
            M.Y(:,k,:) = 1 - bsxfun(@times, (same_context(:,k,:)-diff_context(:,k,:)), 1./baseline);
            denom_factor(:,k,:) = repmat(baseline, 1, size(same_context,1));
        elseif I.divnorm
            assert(~I.normcorr);
            for q = 1:size(diff_context,3)
                xi = same_context(:,k,q)>I.mindenom;
                M.Y(xi,k,q) = diff_context(xi,k,q) ./ same_context(xi,k,q);
                denom_factor(xi,k,q) = same_context(xi,k,q);
                clear xi;
            end
        else
            M.Y(:,k,:) = same_context(:,k,:)-diff_context(:,k,:);
            denom_factor(:,k,:) = 1;
        end
    end
    
    % segment dependent weights
    w_segs = tranweightnsegsfn(L.n_total_segs);
    
    % valid segment durations
    valid_seg_durs = find(w_segs>0);
        
    % integration period
    M.intper_sec = logspace(log10(L.unique_segs(1)/1000), log10(L.unique_segs(end)/1000), 20);
    M.delay_smp = 1:0.5*L.sr;
    M.shape = I.shape;
    M.err = nan(length(M.intper_sec), length(M.delay_smp), length(M.shape), n_channels);
    M.Yh = nan([n_lags, n_seg, length(M.intper_sec), length(M.delay_smp), length(M.shape), n_channels]);
    M.causal_win = false(length(M.intper_sec), length(M.delay_smp), length(M.shape));
    for m = 1:length(M.shape)
        
        fprintf('shape %.2f\n', M.shape(m));
        for i = 1:length(M.intper_sec)
            for j = 1:length(M.delay_smp)
                
                                
                %                 % parameters to run code within just this loop
                %                 I.distr = 'gamma';
                %                 L.lag_t = [0 1];
                %                 L.sr = 100;
                %                 i = 1;
                %                 j = 1;
                %                 m = 1;
                %                 M.intper_sec = 0.1;
                %                 M.delay_smp = L.sr*0;
                %                 M.shape = 1;
                %                 i = 10;
                %                 j = 1;
                %                 m = 1;

                % time vector in samples
                t = -2*L.lag_t(end)*L.sr:2*L.lag_t(end)*L.sr;
                
                switch I.distr
                    case 'gauss'
                        
                        % gaussian window
                        % sigma when int period = central 95%
                        sig_sec = M.intper_sec(i) / 3.92;
                        sig_smp = L.sr * sig_sec;
                        h = normpdf(t - M.delay_smp(j), 0, sig_smp);
                        if I.forcecausal
                            h(t < 0) = 0;
                        end
                        h = h / sum(h);
                        
                        % plot
                        figure;
                        plot(t/L.sr, h);
                        xlim(3*M.intper_sec(i)*[-1, 1] + M.delay_smp(j)/L.sr);
                        
                        % Gaussian window is never causal
                        M.causal_win(i,j,m) = false;
                        
                    case 'gamma'
                        
                        % gaussian window
                        % sigma corresponding to 95%
                        a = M.shape(m);
                        b = 1/M.shape(m);
                        
                        % time vector in samples and seconds
                        t = -2*L.lag_t(end)*L.sr:2*L.lag_t(end)*L.sr;
                        t_sec = t/L.sr;
                        
                        % ratio which to scale stimulus
                        default_intper = gaminv(0.975,a,b) - gaminv(0.025,a,b);
                        r = M.intper_sec(i)/default_intper;
                        
                        % offset to adust delay
                        min_peak = max((a-1)*b,0)*r;
                        c = M.delay_smp(j)/L.sr - min_peak;
                        if M.delay_smp(j)/L.sr < min_peak
                            M.causal_win(i,j,m) = false;
                        else
                            M.causal_win(i,j,m) = true;
                        end
                        
                        % gamma distribution
                        h = gampdf((t_sec-c)/r, a, b);
                        if I.forcecausal
                            h(t_sec < 0) = 0;
                        end
                        h = h / sum(h);
                        
                        % plot
                        % figure;
                        % plot(t_sec, h);
                        % xlim(3*M.intper_sec(i)*[-1, 1] + M.delay_smp/L.sr)
                        
                    otherwise
                        error('No matching distribution');
                end
                
                % calculate predictions for each segment
                for k = valid_seg_durs
                    
                    % rectangular segment
                    seg = zeros(size(t));
                    seg_dur_smp = round(L.sr * L.unique_segs(k)/1000);
                    seg(t <= seg_dur_smp & t >= 0) = 1;
                    % plot(t, [h'/max(h(:)), seg'/max(seg(:))]);
                    
                    % convolve receptive field with segment
                    area = myconv(seg', h', 'causal', false);
                    % plot(t, area);
                    % ylim([0 1]);
                    
                    % prdictions for each electrode
                    for q = 1:n_channels
                        predictor = area(t >= 0);
                        predictor = predictor(1:n_lags);
                        predictor = predictor(:);
                        if I.normcorr || I.divnorm
                            M.Yh(:,k,i,j,m,q) = predictor;
                        else
                            M.Yh(:,k,i,j,m,q) = (1-predictor)*pinv(1-predictor)*M.Y(:,k,q);
                        end
                    end
                end
                
                % calculate error
                for q = 1:n_channels
                    E = (M.Yh(:,valid_seg_durs,i,j,m,q) - M.Y(:,valid_seg_durs,q)).^2;
                    if I.weightdenom
                        w_denom = denom_factor(:,valid_seg_durs,q);
                        w_denom(~(w_denom>I.mindenom)) = 0;
                        w_denom = tranweightdenomfn(w_denom);
                        w_total = bsxfun(@times, w_denom, w_segs(valid_seg_durs));
                    else
                        w_total = w_segs(valid_seg_durs);
                    end
                    w_total = w_total / sum(w_total(:));
                    Ew = bsxfun(@times, E, w_total);
                    clear w_total w_denom;
                    M.err(i,j,m,q) = nanmean(Ew(:));
                end
            end
        end
    end
    
    % find best prediction
    M.Ybest = nan(n_lags, n_seg, n_channels);
    M.best_intper_sec = nan(n_channels,1);
    M.best_delay_smp = nan(n_channels,1);
    M.best_shape_smp = nan(n_channels,1);
    for q = 1:n_channels
        X = M.err(:,:,:,q);
        if any(M.causal_win(:)) && I.bestcausal
            X(~M.causal_win) = inf;
        end
        [~,xi] = min(X(:));
        [a,b,c] = ind2sub(size(X), xi);
        M.Ybest(:,:,q) = M.Yh(:,:,a,b,c,q);
        M.best_intper_sec(q) = M.intper_sec(a);
        M.best_delay_smp(q) = M.delay_smp(b);
        M.best_shape(q) = M.shape(c);
        clear X;
    end
    
    M.channels = L.channels;
    
    save(MAT_file, 'M', '-v7.3');
    
else
    
    load(MAT_file, 'M')
    
end

%% Plot the results

if I.plot_figure
        
    figh = figure;
    
    plot_win_string = [num2str(I.plot_win(1)) '-' num2str(I.plot_win(2))];
    
    n_channels = length(L.channels);
    for q = 1:n_channels
        chan = L.channels(q);
        
        if strcmp(L.chnames{chan}, ['ch' num2str(chan)])
            chname = L.chnames{chan};
        else
            chname = ['ch' num2str(chan) '-' L.chnames{chan}];
        end
        
        % plot prediction for best delay, lineplot
        clf(figh);
        set(figh, 'Position', [100, 100, 900, 900]);
        if I.normcorr || I.divnorm
            corr_range = [-0.5 1.5];
            invariance_line = 1;
        else
            X = M.Y(:,:,q);
            corr_range = quantile(X(:), [0.01, 0.99]);
            clear X;
            invariance_line = NaN;
        end
        valid_seg_durs = find(L.n_total_segs>0);
        for k = valid_seg_durs
            subplot(4, 2, k);
            if I.normcorr
                plot(L.lag_t([1 end]) * 1000, [1,1], 'k--', 'LineWidth', 2);
            end
            hold on;
            plot(L.lag_t([1 end]) * 1000, [0,0], 'k--', 'LineWidth', 2);
            h1 = plot(L.lag_t * 1000, M.Y(:,k,q), 'LineWidth', 2);
            h2 = plot(L.lag_t * 1000, M.Ybest(:,k,q), 'LineWidth', 2);
            plot(L.unique_segs(k)*[1 1], corr_range, 'k--', 'LineWidth', 2);
            if ~isnan(invariance_line); plot(I.plot_win * 1000, invariance_line*[1 1], 'k--', 'LineWidth', 2); end
            xlim(I.plot_win * 1000);
            ylim(corr_range);
            xlabel('Time Lag (ms)');
            ylabel('Pearson Correlation');
            title(sprintf('Seg: %.0f ms', L.unique_segs(k)))
            if k == 1
                legend([h1, h2], 'Data', 'Model');
            end
            box off;
        end
        fname = [L.figure_directory '/' chname '-model-prediction-lineplot_' param_string_modelfit '_plotwin' plot_win_string];
        print_wrapper([fname '.png']);
        print_wrapper([fname '.pdf']);
        
        % plot prediction for best delay, image
        clf(figh);
        set(figh, 'Position', [100, 100, 900, 600]);
        subplot(2,1,1);
        imagesc(M.Y(:,:,q)', corr_range(2) * [-1, 1]);
        subplot(2,1,2);
        imagesc(M.Ybest(:,:,q)', corr_range(2) * [-1, 1]);
        for i = 1:2
            subplot(2,1,i);
            colorbar;
            colormap(flipud(cbrewer('div', 'RdBu', 128)));
            set(gca, 'YTick', 1:length(L.unique_segs), 'YTickLabel', L.unique_segs);
            xticks = get(gca, 'XTick');
            set(gca, 'XTick', xticks, 'XTickLabel', L.lag_t(xticks)*1000);
            ylabel('Seg Dur (ms)'); xlabel('Lag (ms)');
            if i == 2
                title(sprintf('rf=%.f ms, delay=%.f ms', M.best_intper_sec(q)*1000, M.best_delay_smp(q)/L.sr*1000));
            end
        end
        fname = [L.figure_directory '/' chname '-model-prediction_' param_string_modelfit '_plotwin' plot_win_string];
        print_wrapper([fname '.png']);
        % export_fig([fname '.png'], '-png', '-transparent', '-r100');
        
        % plot the error vs. parameters
        X = M.err(:,:,M.best_shape(q)==M.shape,q);
        if any(M.causal_win(:)) && I.plotcausal
            causal_string = '_only-causal';
            X(~M.causal_win(:,:,M.best_shape(q)==M.shape)) = NaN;
        else
            causal_string = '';
        end
        bounds = quantile(-X(:), [0.8, 1]);
        if ~all(isnan(X(:)))
            clf(figh);
            set(figh, 'Position', [100, 100, 600, 600]);
            
            imagesc(-X, bounds);
            colormap(parula(128));
            colorbar;
            ytick = interp1(log2(M.intper_sec), 1:length(M.intper_sec), log2(M.intper_sec));
            set(gca, 'YTick', ytick, 'YTickLabel', M.intper_sec*1000);
            xticks = get(gca, 'XTick');
            set(gca, 'XTickLabel', 1000*M.delay_smp(xticks)/L.sr);
            xlabel('Delay (ms)'); ylabel('Receptive Field Duration (ms)');
            title(sprintf('rf=%.f ms, delay=%.f ms', M.best_intper_sec(q)*1000, M.best_delay_smp(q)/L.sr*1000));
            set(gca, 'FontSize', 12);
            fname = [L.figure_directory '/' chname '-model-error_' param_string_modelfit causal_string];
            print_wrapper([fname '.png']);
            % export_fig([fname '.png'], '-png', '-transparent', '-r100');
            clear X;
        end
    end
    
end