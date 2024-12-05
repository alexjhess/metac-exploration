%% [METAC] exploration data set - task data


%% load data

df = readtable('P:METAC_Iglesias\Data\METAC1_DATA_2024-11-19_1133.csv');

%% extract FAS scores and features
% [df.ppid(find(~isnan(df.fas_exp_total_score))),...
y_fas = df.fas_exp_total_score(find(~isnan(df.fas_exp_total_score)));

X = [df.demo_age(find(~isnan(df.demo_age))),...
    df.demo_gender(find(~isnan(df.demo_gender)))...
    ];

tab = table(df.fas_exp_total_score(find(~isnan(df.fas_exp_total_score))),...
    df.demo_age(find(~isnan(df.demo_age))),...
    df.demo_gender(find(~isnan(df.demo_gender)))...
    );
tab.Properties.VariableNames = {'FAS', 'age', 'gender'};

%% save Y & X
save_path = fullfile('data', ['tmp_data.csv']);
writetable(tab, save_path);


%% find PPIDs
% for n = 1%:size(id,1)
%     ppid = ['TNU_METAC_', num2str(id{n,1})];
%     d = dir(fullfile(basedir, ppid, fname));
%     for i = 1:size(d,1)
%         if contains(d(i).name, 'experiment')
%             % specs stored in rows 1-2
%             dat = readtable(fullfile(d(i).folder, d(i).name), 'NumHeaderLines', 2)
%         end
%     end
% end


