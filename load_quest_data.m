function quest = load_quest_data(dataset)


%% load ppids

if dataset == 1; % 1=discovery set
    id = readtable(fullfile('data', 'metac_ppids_discovery_set.txt'));
elseif dataset == 2 % validation set
    id = readtable(fullfile('data', 'metac_ppids_validation_set.txt'));
end
    
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
y_mfis = df_discovery.mfis_exp_total_score(find(~isnan(df_discovery.mfis_exp_total_score)));
maia3 = df_discovery.maia2_3(find(~isnan(df_discovery.maia2_3)));
maia8 = df_discovery.maia2_8(find(~isnan(df_discovery.maia2_8)));
maia38 = maia3+maia8;
psqi = df_discovery.psqi_score(find(~isnan(df_discovery.psqi_score)));
X = [df_discovery.demo_age(find(~isnan(df_discovery.demo_age))),...
    df_discovery.demo_gender(find(~isnan(df_discovery.demo_gender)))...
    ];

%%
diff_pre_exp = df.exp_diff_perscreening_experiment(find(~isnan(df_discovery.exp_diff_perscreening_experiment)));


%% create data struct
quest = struct;
quest.df_discovery = df_discovery;
quest.y_fas = y_fas;
quest.y_mfis = y_mfis;
quest.maia3 = maia3;
quest.maia8 = maia8;
quest.maia38 = maia38;
quest.psqi = psqi;
quest.X = X;
quest.age = df_discovery.demo_age(find(~isnan(df_discovery.demo_age)));
quest.gender = df_discovery.demo_gender(find(~isnan(df_discovery.demo_gender)));


%% save struct
if dataset == 1 % discovery set
    save('data\discovery_set_quest_tmp.mat', 'quest', '-mat');
elseif dataset == 2 % validation set
    save('data\validation_set_quest_tmp.mat', 'quest', '-mat');
end


end