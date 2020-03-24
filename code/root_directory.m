function y = root_directory

f = filesep;
if exist(['C:' f 'Users' f 'mateo' f 'documents' f ...
          'science' f 'code' f 'integration_quilt'], 'dir') % MLE laptop
    y = ['C:' f 'Users' f 'mateo' f 'documents' f ...
         'science' f 'code' f 'integration_quilt'];
else
    error('No root directory found');

end