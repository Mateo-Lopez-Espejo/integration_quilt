% Concatenate a bunch of stimuli together into one waveform (i.e. for inspecting
% via audacity)

% 2018-03-26: Documented, Sam NH

addpath(genpath([root_directory '/general-analysis-code']));
addpath(genpath([root_directory '/export_fig_v3']));

project_directory = [root_directory '/scrambling-ferrets'];
input_directory = [project_directory '/stimuli/naturalsound-v2'];

stims = mydir(input_directory);

sr = 44100;
dur_per_stim_sec = 0.5; 

stim = nan(sr*dur_per_stim_sec, 1);

for i = 1:length(stims)
    
    [wav, orig_sr] = audioread([input_directory '/' stims{i}]);
    
    if size(wav,2)
        wav = wav(:,1);
    end
    
    if orig_sr ~= sr
        wav = resample(wav, sr, orig_sr);
    end
    
    wav = wav(1:dur_per_stim_sec*sr);
    
    wav = ramp_hann(wav, sr, 0.010);
    
    wav = 0.01*wav/sqrt(mean(wav.^2));
    
    stim((1:dur_per_stim_sec*sr) + (i-1)*dur_per_stim_sec*sr) = wav;
    
end

audiowrite('stim.wav', stim, sr)
% plot(stim)