% CHANGE
% directory with the data
project_directory = [root_directory '/scrambling-ferrets'];
data_directory = [project_directory '/data'];

% CHANGE
% add external repositories
addpath(genpath([root_directory '/general-analysis-code']));
addpath(genpath([root_directory '/export_fig_v3']));

% which recording session to analyze
recording_id = 'tomette012a06_p_NSD';

% smoothing kernel to use
% specified as the FWHM of a Gaussian kernel in units of milliseconds
fwhm_ms = 10;

% whether or not to use spike-sorted data
single_unit = true;

% directory with the recording data
recording_directory = [data_directory '/Tomette/' recording_id(1:10)];

%% Test-retest reliability

% CHANGE
% directory to save figures
figure_directory = [project_directory '/figures/test-retest'];

% test-retest reliability
r = reliability(recording_directory, recording_id, ...
    'plot', true, 'figure_directory', figure_directory, ...
    'fwhm_ms', 10, 'single_unit', single_unit);

%% Lag correlatin analysis

% CHANGE
% directory to save results of analysis
analysis_directory = [project_directory '/analysis/lag-correlation/' recording_id];

% CHANGE
% directory to save figures
figure_directory = [project_directory '/figures/lag-correlation/' recording_id];

% run analysis
% this code is just a wrapper for 
% lag_corr_cross_segdur_modular.m
% and modelfit_lagcorr_cross_segdur_modular.m
[L, M, D, t, S] = lag_corr_cross_segdur(recording_directory, recording_id, ...
    'figure_directory', figure_directory, ...
    'analysis_directory', analysis_directory, ...
    'fwhm_ms', 10, 'single_unit', single_unit);

% save results
save([analysis_directory '/data-and-analysis.mat'], 'L', 'M', 'D', 't', 'S', '-v7.3');