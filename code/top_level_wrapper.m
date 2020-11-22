% Top-level script for all analyses

% 2018-03-18: Created, Sam NH

% paths to external code repositories

root_directory = '/auto/users/mateo/integration_quilt'; % MLE  
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
plot_figure = false;
save_rasters = false; % saves the rasters for Sam, skips the analysis. 

if single_unit
    if lbhb
        recording_ids = {'AMT028b05_p_NTI', 'AMT031a13_p_NTI', 'AMT032a11_p_NTI'...
                         'AMT026a14_p_NTI', 'CRD002a16_p_NTI', 'CRD003b14_p_NTI'...
                         'CRD004a14_p_NTI', 'DRX008b14_p_NTI', 'DRX021a19_p_NTI'...
                         'DRX023a22_p_NTI'};
    else 
        recording_ids = {'tomette002a10_p_NSD'};
    end
    
else
    if lbhb
        recording_ids = {'AMT028b05_p_NTI', 'AMT032a11_p_NTI'};
    else 
        recording_ids = {'tomette002a10_p_NSD'};
    end
end

% asigns full animal name to the animal code 
animal_codes.AMT = 'Amanita';
animal_codes.DRX = 'Drechsler';
animal_codes.CRD = 'Cordyceps';

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

varargin = {};

%% reliability analysis

clc;
if plot_figure
    rall = cell(1, length(recording_ids));
    for i = 1:length(recording_ids)
        try
            fprintf('%d: %s\n', i, recording_ids{i}); drawnow;
            recording_directory = recording_directories{i};
            figure_directory = [project_directory '/figures/test-retest'];
            rall{i} = reliability(recording_directory, recording_ids{i}, ...
                'plot', false, 'figure_directory', figure_directory, ...
                'fwhm_ms', 10, 'single_unit', single_unit, 'spike_thresh', 4);
        catch
            disp([recording_ids{i} ' no such data']);
        end
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
    raster_directory = [project_directory '/rasters/' recording_ids{i}];
    figure_directory = [project_directory '/figures/lag-correlation/' recording_ids{i}];
    r = reliability(recording_directory, recording_ids{i}, ...
        'plot', false, 'figure_directory', figure_directory, ...
        'fwhm_ms', 10, 'single_unit', single_unit, 'spike_thresh', 4);
    lag_corr_cross_segdur(recording_directory, recording_ids{i}, ...
        'figure_directory', figure_directory, ...
        'analysis_directory', analysis_directory, ...
        'raster_directory', raster_directory, ...
        'fwhm_ms', 10, ...
        'single_unit', single_unit, ...
        'units', find(r>0.1),...
        'overwrite', overwrite, ...
        'save_rasters', save_rasters, ...
        'plot_figure', plot_figure,...
        'spike_thresh', 4);
end