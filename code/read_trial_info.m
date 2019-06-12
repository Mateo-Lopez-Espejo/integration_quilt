function T = read_trial_info(recording_directory, recording_id)

% Reads in the stimulus/trial information for a given recording session

% 2018-03-26: Last edited/documented, Sam NH

run([recording_directory '/' recording_id '.m']);

% total numer of trials
T.n_trials = exptevents(end).Trial; %#ok<NODEF>

% exptevents struct into stimulus labels
notes = {exptevents(:).Note};
start_times = cat(1, exptevents(:).StartTime);
trial_indices = cat(1, exptevents(:).Trial);
T.stim_labels = cell(1, T.n_trials);
onset_times = nan(1, T.n_trials);
for i = 1:T.n_trials
    xi = find(trial_indices == i);
    for j = 1:length(xi)
        s = regexp(notes{xi(j)}, 'seg[-\w]*', 'match');
        if ~isempty(s)
            T.stim_labels(i) = s;
        end
        if ~isempty(s) && strcmp(notes{xi(j)}(1:4), 'Stim')
            onset_times(i) = start_times(xi(j));
        end
    end
end

% parse stimulus labels into durations / stims / orders
T.segs = nan(1, T.n_trials);
T.stims = nan(1, T.n_trials);
T.orders = nan(1, T.n_trials);
for i = 1:T.n_trials
    x = regexp(T.stim_labels{i}, '-', 'split');
    T.segs(i) = round(str2double(x{2}(1:end-2)));
    T.stims(i) = round(str2double(x{3}(end)));
    T.orders(i) = round(str2double(x{4}(end)));
end

% come up with a trial order based on the above info
T.trial_order = nan(1, T.n_trials);
T.unique_segs = sort(unique(T.segs));
T.n_seg = length(T.unique_segs);
T.unique_stims = sort(unique(T.stims));
T.n_stims = length(unique(T.stims));
T.unique_orders = sort(unique(T.orders));
T.n_orders = length(unique(T.orders));
count = 0;
for i = 1:T.n_seg
    xi = find(T.unique_segs(i) == T.segs);
    for j = 1:T.n_stims
        yi = find(T.stims(xi) == T.unique_stims(j));
        for k = 1:T.n_orders
            zi = find(T.orders(xi(yi)) == T.unique_orders(k));
            T.trial_order((1:length(zi)) + count) = xi(yi(zi));
            count = count + length(zi);
        end
    end    
end

T.prestim_silence = onset_times(1);
assert(all(abs(onset_times - T.prestim_silence) < 1e-6));