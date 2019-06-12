function y = root_directory

if exist('/home/jennifer', 'dir')
    y = '/home/jennifer/sam';
elseif exist('/home/sam', 'dir')
    y = '/home/sam';
elseif exist('/auto/users/mateo/Sam Analysis', 'dir') % MLE 
    y = '/auto/users/mateo/Sam Analysis';
else
    x = dir('/');
    s = {x(:).name};
    if any(strcmp(s, 'mindhive'))
        y = '/mindhive/nklab/u/svnh';
    else
        if exist('/Users/svnh2/Desktop/projects', 'dir')
            y = '/Users/svnh2/Desktop/projects';
        elseif exist('/Users/scan/Desktop/svnh', 'dir')
            y = '/Users/scan/Desktop/svnh';
        else
            error('No root directory found');
        end
    end
end