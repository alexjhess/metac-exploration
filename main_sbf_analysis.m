%% [METAC] exploration data set - SBF Analysis


%% remove default toolboxes that could be conflicting

% List all paths in search path and break them into a cell array
dirs = regexp(path,['[^;]*'],'match');

% Return all search path directories containing this string
strsToFind = {'tapas', 'hgf-toolbox', 'comb_obs_models', 'VBA-toolbox'}; 

% loop over strings & entries
for i = 1:numel(strsToFind)
    strToFind = char(strsToFind(i));
    % Index to cell entries containing the desired string
    whichCellEntry = find(cellfun(@(dirs) contains(dirs, strToFind), dirs) == 1);
    % remove from path
    for j = 1:length(whichCellEntry)
        rmpath(char(dirs(whichCellEntry(j))))
    end
end

%% add hgf toolbox
addpath(genpath('hgf-toolbox'));
addpath('comb_obs_models');
cd('VBA-toolbox')
VBA_setup();
cd ..

%% load task data

% dat = load_task_discovery_set();
load(fullfile('data', 'discovery_set_tmp.mat'));


%% fit behavioural models & extract other task features

% dat = metac_fit_behaviour(dat);
load(fullfile('data', 'discovery_set_fits_tmp.mat'));

% avg overall control, tolerance, ...
dat.task.avg_c = mean(dat.y_c,'omitnan')';
dat.task.avg_tol = mean(dat.y_tol,'omitnan')';
dat.task.avg_av = mean(dat.y_av,'omitnan')';

%% avg value per sub
figure
plot(dat.task.avg_c)%, '.')
hold on;
plot(dat.task.avg_tol)%, '.')
plot(dat.task.avg_av)%, '.')
hold off;
legend('control', 'tol', 'av')
ylim([0 1])
xlabel('sub')

%% indiv sub
figure
plot(dat.y_c(10:10:80,:))
hold on;
plot(dat.y_tol(10:10:80,:))
plot(dat.y_av(10:10:80,:))
hold off;
ylim([0 1])

%% mean over sub
figure
avg_c = mean(dat.y_c, 2, 'omitnan')
avg_tol = mean(dat.y_tol, 2, 'omitnan')
avg_av = mean(dat.y_av, 2, 'omitnan')
plot(avg_c(10:10:80))
hold on;
plot(avg_tol(10:10:80))
plot(avg_av(10:10:80))
hold off;
ylim([0 1])
xlabel('trial')
legend('control', 'tol', 'av')



%% load quest data

% quest = load_quest_discovery_set();
load(fullfile('data', 'discovery_set_quest_tmp.mat'));

dat.quest = quest;


%% create table for winning mod 1
m = dat.ffx.idx;
pars = dat.param_mat';
if size(pars,2) == 4
    tab = table(dat.quest.y_fas, dat.quest.age, dat.quest.gender, ...
        pars(:,1), pars(:,2), pars(:,3), pars(:,4), ...
        dat.task.avg_c, dat.task.avg_tol, dat.task.avg_av);
    tab.Properties.VariableNames = {'FAS', 'age', 'gender',...
        'gamma', 'shift', 'scale', 'w',...
        'contr', 'tol', 'av'};
end

%% normalize values ?
tab_norm = normalize(tab);


%% write csv file with data for Bayesian ANCOVA
save_path = fullfile('data', ['tmp_data.csv']);
writetable(tab, save_path);

save_path2 = fullfile('data', ['tmp_norm_data.csv']);
writetable(tab_norm, save_path2);

