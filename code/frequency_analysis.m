% Examines the frequency properties of the stimuli and creates whitened
% stimuli.m
% 
% 2018-03-20: Created, Sam NH

addpath(genpath([root_directory '/general-analysis-code']));
addpath(genpath([root_directory '/general-audio-code']));
addpath(genpath([root_directory '/spectrotemporal-synthesis-v2']));

% parameters
P.dur_per_stim_sec = 0.5;
P.audio_sr = 44100;
P.lo_freq_hz = 100;
P.n_filts = round(freq2erb_ferret(P.audio_sr/2)-freq2erb_ferret(P.lo_freq_hz));
P.env_sr = 400;
P.logf_spacing = 1/24;
P.overcomplete = 2;
P.compression_factor = 0.3;
P.animal = 'ferret';
param_idstring = [...
    'dur' num2str(P.dur_per_stim_sec) ...
    '_audsr' num2str(P.audio_sr) ...
    '_lofreq' num2str(P.lo_freq_hz) ...
    '_nfilts' num2str(P.n_filts) ...
    '_envsr' num2str(P.env_sr) ...
    '_freqsr' num2str(1/P.logf_spacing) ...
    '_over' num2str(P.overcomplete) ...
    '_exp' num2str(P.compression_factor) ...
    '_' num2str(P.animal)];

% directory structure
input_directory_name = 'naturalsound20';
project_directory = [root_directory '/scrambling-ferrets'];
input_directory = [project_directory '/stimuli/' input_directory_name];
coch_directory = [project_directory '/analysis/cochleagrams/' ...
    input_directory_name '/' param_idstring];
figure_directory = [project_directory '/figures/cochleagrams/' ...
    input_directory_name '/' param_idstring];
if ~exist(coch_directory, 'dir'); mkdir(coch_directory); end
if ~exist(figure_directory, 'dir'); mkdir(figure_directory); end

% stimuli
stims = mydir(input_directory);

%% Measure excitation patterns

fprintf('Measuring excitation patterns\n\n'); drawnow;

for i = 1:length(stims)
    
    fprintf('%d: %s\n', i, stims{i});
    drawnow;
    
    % read stimulus
    [wav, orig_sr] = audioread([input_directory '/' stims{i}]);
    
    % format
    if size(wav,2)
        wav = wav(:,1);
    end
    if orig_sr ~= P.audio_sr
        wav = resample(wav, P.audio_sr, orig_sr);
    end
    wav = wav(1:P.dur_per_stim_sec*P.audio_sr);
    wav = ramp_hann(wav, P.audio_sr, 0.010);
    wav = 0.01*wav/sqrt(mean(wav.^2));
    
    % cochleogram
    coch_MAT_file = [coch_directory '/' strrep(stims{i}, '.wav', '.mat')];
    if ~exist(coch_MAT_file, 'file')
        [coch, P, R] = wav2coch_without_filts(wav, P);
        save(coch_MAT_file, 'coch', 'P', 'R');
    else
        load(coch_MAT_file,  'coch', 'P');
    end
    
    % store
    if i == 1
        coch_all = nan(size(coch,1), size(coch,2), length(stims));
    end
    coch_all(:,:,i) = coch;
        
end

%% Excitation patterns and excitation standard deviation

addpath(genpath([root_directory '/export_fig_v3']));
excitation_patterns = squeeze_dims(mean(coch_all,1),1);
mean_excitation_pattern = mean(excitation_patterns,2);
figure;
semilogx(P.f, excitation_patterns);
hold on;
semilogx(P.f, mean_excitation_pattern, 'k-', 'LineWidth', 3);
xlim([P.lo_freq_hz, P.audio_sr/2]);
ylim([0 max(excitation_patterns(:))*1.1]);
xlabel('Frequency (Hz)');
ylabel('Compressed Power');
box off;
export_fig([figure_directory '/excitation_patterns.pdf'], '-pdf', '-transparent');

dims = size(coch_all);
coch_unwrap = reshape(permute(coch_all, [1, 3, 2]), [dims(1)*dims(3), dims(2)]);
excitation_std = std(coch_unwrap);
figure;
semilogx(P.f, excitation_std, 'k-', 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Power (std)');
ylim([0 max(excitation_std)*1.1]);
xlim([P.lo_freq_hz, P.audio_sr/2]);
ylim([0 max(excitation_std(:))*1.1]);
box off;
export_fig([figure_directory '/excitation_std.pdf'], '-pdf', '-transparent');

%% Adjust excitation pattern

whitened_stimuli_directory = [project_directory '/stimuli/' input_directory_name '-whitened'];
whitened_coch_directory = [project_directory '/analysis/cochleagrams-whitened/' ...
    input_directory_name '/' param_idstring];
if ~exist(whitened_stimuli_directory, 'dir'); mkdir(whitened_stimuli_directory); end
if ~exist(whitened_coch_directory, 'dir'); mkdir(whitened_coch_directory); end

target_excitation_pattern = ones(size(mean_excitation_pattern)) * mean(mean_excitation_pattern(:));
delta_excitation_pattern = target_excitation_pattern ./ mean_excitation_pattern;

for i = 1:length(stims)
    
    fprintf('%d: %s\n', i, stims{i});
    drawnow;
    
    % create the whitened stimulus
    whitened_stim = [whitened_stimuli_directory '/' stims{i}];
    if ~exist(whitened_stim, 'file')
        
        % load cochleogram
        coch_MAT_file = [coch_directory '/' strrep(stims{i}, '.wav', '.mat')];
        load(coch_MAT_file, 'coch', 'P', 'R');
        
        % whiten
        coch_modified = bsxfun(@times, delta_excitation_pattern', coch);
        
        % convert to waveform
        wav = coch2wav_without_filts(coch_modified, P, R);
        
        % write to file
        audiowrite(whitened_stim, wav, P.audio_sr);
        
    end
    
    % cochleogram
    coch_whitened_MAT_file = [whitened_coch_directory '/' strrep(stims{i}, '.wav', '.mat')];
    if ~exist(coch_whitened_MAT_file, 'file')
        wav = audioread(whitened_stim);
        [coch, P, R] = wav2coch_without_filts(wav, P);
        save(coch_whitened_MAT_file, 'coch', 'P', 'R');
    else
        load(coch_whitened_MAT_file,  'coch');
    end
    
    % store
    if i == 1
        coch_whitened_all = nan(size(coch,1), size(coch,2), length(stims));
    end
    coch_whitened_all(:,:,i) = coch;
    
end

%% Excitation patterns and standard deviation for whitened stimuli

addpath(genpath([root_directory '/export_fig_v3']));
whitened_excitation_patterns = squeeze_dims(mean(coch_whitened_all,1),1);
whitened_mean_excitation_pattern = mean(whitened_excitation_patterns,2);
figure;
semilogx(P.f, whitened_excitation_patterns);
hold on;
semilogx(P.f, whitened_mean_excitation_pattern, 'k-', 'LineWidth', 3);
xlim([P.lo_freq_hz, P.audio_sr/2]);
ylim([0 max(whitened_excitation_patterns(:))*1.1]);
xlabel('Frequency (Hz)');
ylabel('Compressed Power');
box off;
export_fig([figure_directory '/excitation_patterns_whitened.pdf'], ...
    '-pdf', '-transparent');

dims = size(coch_whitened_all);
coch_whitened_unwrap = reshape(permute(coch_whitened_all, [1, 3, 2]), [dims(1)*dims(3), dims(2)]);
excitation_whitened_std = std(coch_whitened_unwrap);
figure;
semilogx(P.f, excitation_whitened_std, 'k-', 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Power (std)');
ylim([0 max(excitation_whitened_std)*1.1]);
xlim([P.lo_freq_hz, P.audio_sr/2]);
ylim([0 max(excitation_whitened_std(:))*1.1]);
box off;
export_fig([figure_directory '/excitation_std_whitened.pdf'], ...
    '-pdf', '-transparent');