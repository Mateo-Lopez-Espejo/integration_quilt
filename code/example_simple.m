% load data and analyses for one recording

recording_id = 'tomette012a06_p_NSD';
analysis_directory = [root_directory '/scrambling-ferrets/analysis/lag-correlation/' recording_id];
load([analysis_directory '/data-and-analysis.mat'], 'L', 'M', 'D', 't', 'S');

%% lag correlation analysis

L_new = lag_corr_cross_segdur_modular(D, t, S, 'overwrite', true);

%% Fit lag correlations with model

M_new = modelfit_lagcorr_cross_segdur_modular(L_new, 'overwrite', true, 'divnorm', true, 'weightdenom', true);