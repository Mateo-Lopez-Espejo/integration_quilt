function y = data_directory

f = filesep;
if exist(['H:' f 'daq'], 'dir') % MLE laptop 
    y = ['H:' f 'daq'];
else
    error('No data directory found');
end