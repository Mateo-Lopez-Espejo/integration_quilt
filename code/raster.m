function [raster_smooth, tbins] = raster(recording_directory, recording_id, win, varargin)

% Computes a raster for a given recording. See first cell for key optional
% parameters.

% 2018-03-26: Last edited

%% Parameters

I.fwhm_ms = 10;
I.raster_sr = 500;
I.spike_thresh = -4;
I.spike_sr = 31250;
I.single_unit = false;
I = parse_optInputs_keyvalue(varargin, I);

% determine number of trials
if I.single_unit
    n_trials = load([recording_directory '/sorted/' recording_id '.spk.mat'],...
                    'nrec');
    n_trials = n_trials.nrec;
else
    elec = 1;
    load([recording_directory '/tmp/' recording_id '.001.1.elec' ...
        num2str(elec) '.sig' num2str(I.spike_thresh) '.NOCOM.mat'], ...
        'trialid', 'spikebin');
    n_trials = max(trialid);
end

% count number of units
tbins = win(1) : 1/I.raster_sr : win(2);
n_bins = length(tbins);

% raster from multi or single unit
if I.single_unit
    
    % raster for all units
    raster = nan(n_bins, n_trials, 100);
    X = load([recording_directory '/sorted/' recording_id '.spk.mat']);
    if isfield(X, 'rate') && ~isempty(X.rate)% MLE. Whis is this value coming from outside?
        I.spike_sr = X.rate; 
    end
    n_electrodes = length(X.sortinfo);
    n_units = 0;
    for elec = 1:n_electrodes
        if ~isempty(X.sortinfo{elec})
            for j = 1:length(X.sortinfo{elec}{1})
                % lbhb compatibility with sorted spikes data structures.
                % skips over empty units. MLE 2019 06 11
                if ~isempty(X.sortinfo{elec}{1}(j).unitSpikes)
                    spike_data = X.sortinfo{elec}{1}(j).unitSpikes;
                    n_units = n_units + 1;
                    for i = 1:n_trials
                        spike_times = spike_data(2, spike_data(1, :)==i) / I.spike_sr;
                        raster(:, i, n_units) = myhist(spike_times, tbins);
                    end
                end
            end
        end
    end
    raster = raster(:,:,1:n_units);
        
else
    
    % determine number of electrodes
    elec = 0;
    while 1
        elec = elec+1;
        fname = [recording_directory '/tmp/' recording_id '.001.1.elec' ...
            num2str(elec) '.sig' num2str(I.spike_thresh) '.NOCOM.mat'];
        if ~exist(fname,'file')
            n_electrodes = elec-1;
            break;
        end
    end
    
    % bins for the raster
    tbins = win(1) : 1/I.raster_sr : win(2);
    n_bins = length(tbins);
    
    % raster for all units
    raster = nan(n_bins, n_trials, n_electrodes);
    for elec = 1:n_electrodes
        load([recording_directory '/tmp/' recording_id '.001.1.elec' ...
            num2str(elec) '.sig' num2str(I.spike_thresh) '.NOCOM.mat'], ...
            'trialid', 'spikebin');
        for i = 1:n_trials
            spike_times = spikebin(trialid == i)/I.spike_sr;
            if ~isempty(spike_times)
                raster(:, i, elec) = myhist(spike_times, tbins);
            end
        end
    end
end

% smooth with Gaussian kernel
fwhm_ms = 10;
if fwhm_ms > 0
    raster_smooth = mysmooth(raster, I.raster_sr, I.fwhm_ms);
else
    raster_smooth = raster;
end