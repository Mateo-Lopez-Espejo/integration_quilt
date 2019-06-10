clc;

% paths to external code repositories
% addpath(genpath([root_directory '/general-analysis-code']));
% addpath(genpath([root_directory '/export_fig_v3']));
addpath(genpath([root_directory '/baphy']));
rmpath([root_directory '/baphy/Utilities/CommonCodes']);
baphy_set_path;

%%

clc;
baphy_remote;

%%

mfilenamc = '/Volumes/data/Tomette/tomette014/tomette014a08_p_TOR.m';
load(mfilenamc)

%%

MD_computeCSD('Identifier', 'tomette014a08_p_TOR', 'Electrodes', [1:32]);


% project_directory = [root_directory '/scrambling-ferrets'];




figure;
mfilenamc = '/Volumes/data/Tomette/tomette014/tomette014a08_p_TOR.m';
spikefile = '/Volumes/data/Tomette/tomette014/sorted/tomette014a08_p_TOR.spk.mat';
channel = 1;
unit = 1;
axeshandle = gca;
STRFestimation(mfilenamc,spikefile,channel,unit,axeshandle)
