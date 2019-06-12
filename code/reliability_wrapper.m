% Illustrates how to call reliability.m

% 2018-03-26: Last edited

% recording_id = 'tomette014a09_p_NSD';
% recording_directory = '/Volumes/data/Tomette/tomette014';
recording_id = 'tomette012a06_p_NSD';
recording_directory = '/Volumes/data/Tomette/tomette012';
% recording_id = 'tomette002a10_p_NSD';
% recording_directory = '/Volumes/data/Tomette/tomette002';


recording_id = 'tomette012a06_p_NSD';
recording_directory = '/Volumes/data/Tomette/tomette012';
reliability(recording_directory, recording_id, 'n_reps', 4, 'individ_electrodes', [4, 20])