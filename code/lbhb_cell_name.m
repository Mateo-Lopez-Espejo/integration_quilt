function [cellnames] = lbhb_cell_name(recording_directory, recording_id)
%LBHB_CELL_NAME given a .spk.mat file returns the singe unit names 
%following LBHB format containing site electrode and cell e.g. 'AMT032a-14-1'
%the order corresponds to the oder of single unis as given by Sam Raster
%function
% 2019-06-11: Created, MLE
%   Detailed explanation goes here

X = load([recording_directory '/sorted/' recording_id '.spk.mat']);

site_name = recording_id(1:7);
n_electrodes = length(X.sortinfo);
n_units = 0;
cellnames = cell(1,100);
for elec = 1:n_electrodes
    if ~isempty(X.sortinfo{elec})
        for j = 1:length(X.sortinfo{elec}{1})
            if ~isempty(X.sortinfo{elec}{1}(j).unitSpikes)
                n_units = n_units + 1;
                formatSpec = '%s-%02d-%d';
                cellnames{n_units} = sprintf(formatSpec, site_name, elec, j);
            end
        end
    end
end

cellnames = cellnames(1:n_units);
end

