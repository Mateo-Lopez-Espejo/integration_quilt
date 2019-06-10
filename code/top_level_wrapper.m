% Top-level script for all analyses

% 2018-03-18: Created, Sam NH

% paths to external code repositories

root_directory = '/home/mateo/Sam Analysis'; % MLE  
addpath(genpath([root_directory '/general-analysis-code']));
addpath(genpath([root_directory '/export_fig_v3']));
project_directory = [root_directory '/scrambling-ferrets'];

if exist([project_directory '/data'], 'dir')
    data_directory = [project_directory '/data'];
elseif exist('/auto/data', 'dir')
    data_directory = root_directory; % MLE
    %data_directory = '/auto/data'; % SNH
else
    data_directory = '/Volumes/data';
end

single_unit = true;
mateo = false;

if single_unit
    if mateo
        recording_ids = {'AMT032a11_p_NTI'}; % MLE LBHB example
    else 
        recording_ids = {'tomette002a10_p_NSD'};
    end
    %recording_ids = {...
    %    'tomette002a10_p_NSD', 'tomette012a06_p_NSD', 'tomette014a09_p_NSD', ...
    %    'tomette016b09_p_NSD', 'tomette018a06_p_NSD', ...
    %    'tomette032a07_p_NSD', 'tomette038a06_p_NSD', ...
    %    'tomette063a08_p_NSD', ...
    %    }; % SHN
else
    recording_ids = {'tomette067a17_p_NSD'};
    %recording_ids = {...
    %    'tomette002a10_p_NSD', 'tomette012a06_p_NSD', 'tomette014a09_p_NSD', ...
    %    'tomette016b09_p_NSD', 'tomette018a06_p_NSD', 'tomette026a08_p_NSD', ...
    %    'tomette032a07_p_NSD', 'tomette038a06_p_NSD', 'tomette048a03_p_NSD', ...
    %    'tomette063a08_p_NSD', 'tomette067a17_p_NSD', ...
    %    }; % SHN
end


recording_directories = cell(1, length(recording_ids));

for i = 1:length(recording_ids)
    if mateo
        recording_directories{i} = [data_directory '/Amanita/' recording_ids{i}(1:6)];
    else
        recording_directories{i} = [data_directory '/Tomette/' recording_ids{i}(1:10)];
    end
end

% recording_directory = [data_directory '/Tomette/tomette002'];
% recording_id = 'tomette002a10_p_NSD';
% recording_directory = [data_directory '/Tomette/tomette012'];
% recording_id = 'tomette012a06_p_NSD';
% recording_directory = [data_directory '/Tomette/tomette014'];
% recording_id = 'tomette014a09_p_NSD';
% recording_directory = [data_directory '/Tomette/tomette016'];
% recording_id = 'tomette016b09_p_NSD';
% recording_directory = [data_directory '/Tomette/tomette018'];
% recording_id = 'tomette018a06_p_NSD';
% recording_directory = [data_directory '/Tomette/tomette032'];
% recording_id = 'tomette032a07_p_NSD';
% recording_directory = [data_directory '/Tomette/tomette038'];
% recording_id = 'tomette038a06_p_NSD';
% recording_directory = [data_directory '/Tomette/tomette048'];
% recording_id = 'tomette048a03_p_NSD';
% recording_directory = [data_directory '/Tomette/tomette063'];
% recording_id = 'tomette063a08_p_NSD';
% recording_directory = [data_directory '/Tomette/tomette067'];
% recording_id = 'tomette067a17_p_NSD';
% recording_directory = [data_directory '/Tomette/tomette026'];
% recording_id = 'tomette026a08_p_NSD';

varargin = {};

%% reliability analysis

clc;
rall = cell(1, length(recording_ids));
for i = 1:length(recording_ids)
    fprintf('%d: %s\n', i, recording_ids{i}); drawnow;
    recording_directory = recording_directories{i};
    figure_directory = [project_directory '/figures/test-retest'];
    rall{i} = reliability(recording_directory, recording_ids{i}, ...
        'plot', false, 'figure_directory', figure_directory, ...
        'fwhm_ms', 10, 'single_unit', single_unit);
end

%% Histogram of test-retest reliability

rvec = cat(2, rall{:});
hist(rvec);
xlabel('test-retest reliability');
ylabel('number of single units');
%export_fig([figure_directory '/hist-allunits-singleunit' num2str(single_unit) '.pdf']);

%% Lag analysis

clc;
for i = 1:length(recording_ids)
    recording_directory = recording_directories{i};
    analysis_directory = [project_directory '/analysis/lag-correlation/' recording_ids{i}];
    figure_directory = [project_directory '/figures/lag-correlation/' recording_ids{i}];
    r = reliability(recording_directory, recording_ids{i}, ...
        'plot', false, 'figure_directory', figure_directory, ...
        'fwhm_ms', 10, 'single_unit', single_unit);
    lag_corr_cross_segdur(recording_directory, recording_ids{i}, ...
        'figure_directory', figure_directory, ...
        'analysis_directory', analysis_directory, ...
        'fwhm_ms', 10, 'single_unit', single_unit, ...
        'units', find(r>0.1));
end