% Top-level script for all analyses

% 2018-03-18: Created, Sam NH

% paths to external code repositories

root_directory = '/auto/users/mateo/Sam Analysis'; % MLE  
addpath(genpath([root_directory '/general-analysis-code']));
addpath(genpath([root_directory '/export_fig_v3']));
addpath(genpath([root_directory '/helper-functions']));
project_directory = [root_directory '/scrambling-ferrets'];

% lbhb file reading utilities, hopefully it does not brek pathing
baphy_set_path

if exist([project_directory '/data'], 'dir')
    %data_directory = [project_directory '/data']; % SNH
elseif exist('/auto/data/daq', 'dir')
    data_directory = '/auto/data/daq';
    %data_directory = '/auto/data'; % SNH
else
    %data_directory = '/Volumes/data'; % SNH
end

single_unit = true;
lbhb = true; % uses data from david lab.
overwrite = false;
plot_figure = true;

if single_unit
    if lbhb
        recording_ids = {'AMT026a14_p_NTI', 'AMT028b05_p_NTI', 'AMT032a11_p_NTI'};
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
    if lbhb
        recording_ids = {'AMT028b05_p_NTI', 'AMT032a11_p_NTI'};
    else 
        recording_ids = {'tomette002a10_p_NSD'};
    end
    %recording_ids = {...
    %    'tomette002a10_p_NSD', 'tomette012a06_p_NSD', 'tomette014a09_p_NSD', ...
    %    'tomette016b09_p_NSD', 'tomette018a06_p_NSD', 'tomette026a08_p_NSD', ...
    %    'tomette032a07_p_NSD', 'tomette038a06_p_NSD', 'tomette048a03_p_NSD', ...
    %    'tomette063a08_p_NSD', 'tomette067a17_p_NSD', ...
    %    }; % SHN
end

% asigns full animal name to the animal code 
animal_codes.AMT = 'Amanita';

recording_directories = cell(1, length(recording_ids));

for i = 1:length(recording_ids)
    a_code = recording_ids{i}(1:3);
    if lbhb
        recording_directories{i} = [data_directory '/'...
                                    animal_codes.(a_code) '/' ...
                                    recording_ids{i}(1:6)];
    else
        recording_directories{i} = ['/auto/users/mateo/Sam Analysis/Tomette/' recording_ids{i}(1:10)];
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
if plot_figure
    rall = cell(1, length(recording_ids));
    for i = 1:length(recording_ids)
        fprintf('%d: %s\n', i, recording_ids{i}); drawnow;
        recording_directory = recording_directories{i};
        figure_directory = [project_directory '/figures/test-retest'];
        rall{i} = reliability(recording_directory, recording_ids{i}, ...
            'plot', false, 'figure_directory', figure_directory, ...
            'fwhm_ms', 10, 'single_unit', single_unit, 'spike_thresh', 4);
    end
end

%% Histogram of test-retest reliability

if plot_figure
    rvec = cat(2, rall{:});
    hist(rvec);
    xlabel('test-retest reliability');
    ylabel('number of single units');
    export_fig([figure_directory '/hist-allunits-singleunit' num2str(single_unit) '.pdf']);
end
    
%% Lag analysis

clc;
for i = 1:length(recording_ids)
    recording_directory = recording_directories{i};
    analysis_directory = [project_directory '/analysis/lag-correlation/' recording_ids{i}];
    figure_directory = [project_directory '/figures/lag-correlation/' recording_ids{i}];
    r = reliability(recording_directory, recording_ids{i}, ...
        'plot', false, 'figure_directory', figure_directory, ...
        'fwhm_ms', 10, 'single_unit', single_unit, 'spike_thresh', 4);
    lag_corr_cross_segdur(recording_directory, recording_ids{i}, ...
        'figure_directory', figure_directory, ...
        'analysis_directory', analysis_directory, ...
        'fwhm_ms', 10, ...
        'single_unit', single_unit, ...
        'units', find(r>0.1),...
        'overwrite', overwrite, ...
        'plot_figure', plot_figure,...
        'spike_thresh', 4);
end