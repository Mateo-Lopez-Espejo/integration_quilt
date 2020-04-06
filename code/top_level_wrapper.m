% Top-level script for all analyses

% 2018-03-18: Created, Sam NH

% paths to external code repositories
f = filesep;

disp(['root directoy: ' root_directory()])
disp( ['data directory: ' data_directory])

addpath(genpath([root_directory() f 'general-analysis-code']));
addpath(genpath([root_directory() f 'export_fig_v3']));
addpath(genpath([root_directory() f 'helper-functions']));
project_directory = [root_directory f 'scrambling-ferrets'];


single_unit = true;
lbhb = true; % uses data from david lab.
overwrite = false;
plot_figure = true;
save_rasters = false; % saves the rasters for Sam, skips the analysis. 

if single_unit
    if lbhb
        recording_ids = {'AMT026a14_p_NTI', 'AMT028b05_p_NTI',...
            'AMT031a13_p_NTI', 'AMT032a11_p_NTI', ...
            'DRX008b14_p_NTI', 'DRX021a19_p_NTI'};
        %recording_ids = {'DRX021a19_p_NTI'};
        %recording_ids = {'AMT026a14_p_NTI'};
        % to sort: 'AMT029a16_p_NTI', AMT030a12_p_NTI, ,
        % 'DRX010c05_p_NTI', 'DRX023a22_p_NTI'
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

recording_directories = cell(1, length(recording_ids));

for i = 1:length(recording_ids)
    a_code = recording_ids{i}(1:3);
    if lbhb
        recording_directories{i} = [data_directory f...
                                    animal_codes.(a_code) f ...
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
            figure_directory = [project_directory f 'figures' f 'test-retest'];
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
    filename = [figure_directory f 'hist-allunits-singleunit' num2str(single_unit) '.pdf'];
    fname = [figure_directory f 'hist-allunits-singleunit' num2str(single_unit)];
    %export_fig(filename);
    print_wrapper([fname '.png']);
end
    
%% Lag analysis

clc;
for i = 1:length(recording_ids)
    disp(['working on: ' recording_ids{i}]);
    recording_directory = recording_directories{i};
    analysis_directory = [project_directory f 'analysis' f ...
                          'lag-correlation' f recording_ids{i}];
    raster_directory = [project_directory f 'rasters' f recording_ids{i}];
    figure_directory = [project_directory f 'figures' f ...
                        'lag-correlation' f recording_ids{i}];
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