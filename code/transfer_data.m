% paths to external code repositories
addpath(genpath([root_directory '/general-analysis-code']));
addpath(genpath([root_directory '/export_fig_v3']));
project_directory = [root_directory '/scrambling-ferrets'];

% directories
transfer_directory = [project_directory '/data/transfer'];
local_directory = [project_directory '/data'];
if ~exist(local_directory, 'dir'); mkdir(local_directory); end
if ~exist(transfer_directory, 'dir'); mkdir(transfer_directory); end

biguy_address = 'sam@129.199.80.246';
biguy_directory = '/auto/data';
kerfuffle_address = 'svnh@kerfuffle.mit.edu';
kerfuffle_directory = '/mindhive/nklab/u/svnh/scrambling-ferrets/data';

files = {...
    'Tomette/tomette014/sorted/tomette014a09_p_NSD.spk.mat', ...
    'Tomette/tomette012/sorted/tomette012a06_p_NSD.spk.mat', ...
    'Tomette/tomette002/sorted/tomette002a10_p_NSD.spk.mat', ...
    'Tomette/tomette016/sorted/tomette016b09_p_NSD.spk.mat', ...
    'Tomette/tomette018/sorted/tomette018a06_p_NSD.spk.mat', ...
    'Tomette/tomette032/sorted/tomette032a07_p_NSD.spk.mat', ...
    'Tomette/tomette038/sorted/tomette038a06_p_NSD.spk.mat', ...
    'Tomette/tomette026/sorted/tomette026a08_p_NSD.spk.mat', ...
    'Tomette/tomette063/sorted/tomette063a08_p_NSD.spk.mat', ...
    'Tomette/tomette002/tomette002a10_p_NSD.m', ...
    'Tomette/tomette012/tomette012a06_p_NSD.m', ...
    'Tomette/tomette014/tomette014a09_p_NSD.m', ...
    'Tomette/tomette016/tomette016b09_p_NSD.m', ...
    'Tomette/tomette018/tomette018a06_p_NSD.m', ...
    'Tomette/tomette032/tomette032a07_p_NSD.m', ...
    'Tomette/tomette038/tomette038a06_p_NSD.m', ...
    'Tomette/tomette048/tomette048a03_p_NSD.m', ...
    'Tomette/tomette063/tomette063a08_p_NSD.m', ...
    'Tomette/tomette067/tomette067a17_p_NSD.m', ...
    'Tomette/tomette026/tomette026a08_p_NSD.m', ...
    };

dirs = {...
    'Tomette/tomette002/tmp', ...
    'Tomette/tomette012/tmp', ...
    'Tomette/tomette014/tmp', ...
    'Tomette/tomette016/tmp', ...
    'Tomette/tomette018/tmp', ...
    'Tomette/tomette026/tmp', ...
    'Tomette/tomette032/tmp', ...
    'Tomette/tomette038/tmp', ...
    'Tomette/tomette048/tmp', ...
    'Tomette/tomette063/tmp', ...
    'Tomette/tomette067/tmp', ...
    };

%% Transfer files to local

clc;
for i = 1:length(files)
    fprintf('%d: %s\n', i, files{i}); drawnow;
    cmd = ['/usr/local/bin/sshpass -p ''matosM@M#'' rsync -r "' ...
        biguy_address ':' biguy_directory '/' files{i} '" ' ...
        mkpdir([local_directory '/' files{i}])];
    unix(cmd);
end

%% Transfer files to kerfuffle

clc;
for i = 1:length(files)
    [p,f,e] = fileparts(files{i});
    fprintf('%d: %s\n', i, files{i}); drawnow;
    cmd = [...
        '/usr/local/bin/sshpass -p "fjdkM&M*" ' ...
        'ssh -o StrictHostKeyChecking=no ' kerfuffle_address ...
        ' mkdir -p ' kerfuffle_directory '/' p];
    unix(cmd);
    cmd = [...
        '/usr/local/bin/sshpass -p "fjdkM&M*" rsync ' ...
        '-e "ssh -o StrictHostKeyChecking=no" ' ...
        local_directory '/' files{i} ' ' ...
        kerfuffle_address ':' kerfuffle_directory '/' files{i}];
    unix(cmd);
end

%% Transfer directories to local

clc;
for i = 6%1:length(dirs)
    fprintf('%d: %s\n', i, dirs{i}); drawnow;
    [pdir, ~, ~] = fileparts(dirs{i});
    biguy_address = 'sam@129.199.80.246';
    cmd = ['/usr/local/bin/sshpass -p ''matosM@M#'' rsync -r "' ...
        biguy_address ':' biguy_directory '/' dirs{i} '" ' ...
        mkpdir([local_directory '/' pdir])];
    unix(cmd);
end

%% Remove files that are not needed

for i = 1:length(dirs)
    all_files = mydir([local_directory '/' dirs{i}]);
    nsd_files = mydir([local_directory '/' dirs{i}], 'NSD');
    files_to_delete = setdiff(all_files, nsd_files);
    for j = 1:length(files_to_delete)
        delete([local_directory '/' dirs{i} '/' files_to_delete{j}]);
    end
end

%% Transfer directories to kerfuffle

clc;
for i = 1:length(dirs)
    fprintf('%d: %s\n', i, dirs{i}); drawnow;
    [pdir, ~, ~] = fileparts(dirs{i});
    cmd = [...
        '/usr/local/bin/sshpass -p "fjdkM&M*" ' ...
        'ssh -o StrictHostKeyChecking=no ' kerfuffle_address ...
        ' mkdir -p ' kerfuffle_directory '/' pdir];
    unix(cmd);
    cmd = ['/usr/local/bin/sshpass -p "fjdkM&M*" rsync ' ...
        '-e "ssh -o StrictHostKeyChecking=no" -r ' ...
        local_directory '/' dirs{i} ' ' ...
        kerfuffle_address ':' kerfuffle_directory '/' pdir];
    unix(cmd);
end


%% remove files that are not needed

% for i = 1%:length(fnames_tmp)
%     mydir([transfer_directory ]);
% end