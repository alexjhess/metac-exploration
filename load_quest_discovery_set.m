function quest = load_quest_discovery_set()


%% load ppids

id = readtable(fullfile('data', 'metac_ppids_exploration_set.txt'));
ppidstr = cell(size(id,1),1);
for i = 1: size(id,1)
    ppidstr{i} = ['METAC_', num2str(id{i,1})];
end


%% load questionnaire data

datadir = fullfile('P:METAC_Iglesias', 'Data');
% datadir = fullfile('data');
f = dir(fullfile(datadir, 'METAC1*.csv'));

%% load ppids of discovery set
df = readtable(fullfile(f.folder, f.name));
idx = find(contains(df.ppid, ppidstr));
df_discovery = df(idx,:);

%% extract FAS scores and features
y_fas = df_discovery.fas_exp_total_score(find(~isnan(df_discovery.fas_exp_total_score)));
X = [df_discovery.demo_age(find(~isnan(df_discovery.demo_age))),...
    df_discovery.demo_gender(find(~isnan(df_discovery.demo_gender)))...
    ];



%% create data struct
quest = struct;
quest.df_discovery = df_discovery;
quest.y_fas = y_fas;
quest.X = X;
quest.age = df_discovery.demo_age(find(~isnan(df_discovery.demo_age)));
quest.gender = df_discovery.demo_gender(find(~isnan(df_discovery.demo_gender)));


%% save struct
save('data\discovery_set_quest_tmp.mat', 'quest', '-mat');


end