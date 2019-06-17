% script to run and plot fited values for all cells analyzes 

% 2019-06-16: Created. Mateo Lopez Espejo

% path to data

root_directory = '/auto/users/mateo/Sam Analysis'; % MLE  
addpath(genpath([root_directory '/general-analysis-code']));
addpath(genpath([root_directory '/export_fig_v3']));
results_dir = [root_directory '/scrambling-ferrets/analysis/lag-correlation'];

recording_ids = {'AMT028b05_p_NTI', 'AMT032a11_p_NTI',  'tomette002a10_p_NSD'};

% finds and lists the directory of and site of relevant sites
all_analysis_files = subdir(results_dir);
files_to_load = cell(1, length(all_analysis_files));
file_counter = 0;
for ff = 1:length(all_analysis_files)
    [filepath, name, ext] = fileparts(all_analysis_files(ff).name);
    for rr=1:length(recording_ids)
        if ~isempty(strfind(filepath, recording_ids{rr}))
            if startsWith(name, 'model_fit_')
                file_counter = file_counter + 1;
                files_to_load{file_counter} = all_analysis_files(ff).name;
            end
        end
    end
end
files_to_load = files_to_load(1:file_counter);


% loads and parses relevant data

best_fits(200).best_intper_sec = [];
best_fits(200).best_delay_smp = [];
best_fits(200).chname = [];
unit_counter = 0 ;

for ff=1:length(files_to_load)
    load(files_to_load{ff}, 'M');
    
    % iterates over units in site
    for uu = 1:length(M.channels)
        unit_counter = unit_counter + 1;
        best_fits(unit_counter).best_intper_sec = M.best_delay_smp(uu)/M.sr*1000;
        best_fits(unit_counter).best_delay_smp = M.best_intper_sec(uu)*1000;
        if isfield(M, 'chnames')
            best_fits(unit_counter).chnames = M.chnames{M.channels(uu)};
        else
            parts = split(files_to_load{ff}, '/');
            best_fits(unit_counter).chnames = [parts{9} '_ch_' M.channels(uu)];
        end
    end
end
best_fits = best_fits(1:unit_counter);

int = zeros(1, length(best_fits));
lag = zeros(1, length(best_fits));
cellname = cell(1, length(best_fits));
for cc = 1:length(best_fits)
    int(cc) = best_fits(cc).best_intper_sec;
    lag(cc) = best_fits(cc).best_delay_smp;
    cellname{cc} = best_fits(cc).chnames;
end 

fig = figure();
scatter(int, lag, 'filled');
dx = 0; dy = 0; % displacement so the text does not overlay the data points
snapnow;
text(int+dx, lag+dy, cellname, 'fontsize', 15);
xlabel('best time Lag (ms)');
ylabel('best integration (ms)');
ylim([0, 510]);
set(gca, 'YDir','reverse')



% parses into convenient array/structure
% plots scatter