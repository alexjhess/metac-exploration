function dat = load_task_validation_set()


%% load behavioural data

basedir = fullfile('R:METAC', 'behavior', 'raw');
fname = fullfile('behavior', 'task', '*.csv');

id = readtable(fullfile('data', 'metac_ppids_exploration_set.txt'));

u_bin = NaN(80,size(id,1));
y_pred = NaN(80,size(id,1));
y_mc = NaN(80,size(id,1));
y_c = NaN(80,size(id,1));
y_tol = NaN(80,size(id,1));
y_av = NaN(80,size(id,1));

for n = 1:size(id,1)
    ppid = ['TNU_METAC_', num2str(id{n,1})];
    disp(ppid)
    d = dir(fullfile(basedir, ppid, fname));
    for i = 1:size(d,1)
        if contains(d(i).name, 'experiment')
            % specs stored in rows 1-2
            tmp = readtable(fullfile(d(i).folder, d(i).name), 'NumHeaderLines', 2);
            if ismember('overall_control', tmp.Properties.VariableNames) %task v4
                tmp_av = readtable(fullfile(d(i).folder, d(i).name), 'NumHeaderLines', 83);
            end
        end
    end
    u_bin(:,n) = tmp.jSuccess(1:80);
    y_pred(:,n) = tmp.prediction(1:80);
    % try
    %     y_mc(:,n) = tmp.control(1:80);
    % catch
    %     y_mc(:,n) = tmp.expl_control(1:80);
    % end
    if ismember('overall_control', tmp.Properties.VariableNames) %task v4
        y_mc(:,n) = tmp.expl_control(1:80);
        y_c(tmp.overall_control(1:80)~=-1,n) = tmp.overall_control(tmp.overall_control(1:80)~=-1);
        y_av(80,n) = tmp_av.aversiveness;
    else % task v3.1
        y_mc(:,n) = tmp.control(1:80);
        y_c(tmp.handling(1:80)~=-1,n) = tmp.handling(tmp.handling(1:80)~=-1);
        y_av(tmp.aversiveness(1:80)~=-1,n) = tmp.aversiveness(tmp.aversiveness(1:80)~=-1);
    end
    y_tol(tmp.tolerance(1:80)~=-1,n) = tmp.tolerance(tmp.tolerance(1:80)~=-1);
end

u_pe = u_bin - y_pred;
u_pe_sq = u_pe.^2;
u_pe_pos = NaN(size(u_pe));
u_pe_pos(u_pe>=0) = u_pe(u_pe>=0);
u_pe_neg = NaN(size(u_pe));
u_pe_neg(u_pe<0) = abs(u_pe(u_pe<0));



%% load experiment structure
exp_seq = readtable(fullfile('data', 'metac_experiment_seq.txt'));

% mab input sequence (wind / no wind)
u_mab2 = NaN(size(exp_seq,1), 1);
u_mab4 = NaN(size(exp_seq,1), 1);

% phases
trials = 1:size(u_bin,1);
phase = ones(size(trials)); % 1: easy phase

n=1;
% different input seq for wind / no wind
for k = 1:size(exp_seq,1)
    if exp_seq.wind(k)==0
        u_mab2(k,1) = 1; % no wind
        if exp_seq.iWidth(k)==0 % small island
            u_mab4(k,1) = 2;
        elseif exp_seq.iWidth(k)==1 % large island
            u_mab4(k,1) = 1;
        end
    elseif exp_seq.wind(k)==1 || exp_seq.wind(k)==-1
        u_mab2(k,1) = 2; % wind
        if exp_seq.iWidth(k)==0 % small island
            u_mab4(k,1) = 4;
        elseif exp_seq.iWidth(k)==1 % large island
            u_mab4(k,1) = 3;
        end
    end
    if k > 10 && k <= 30
        phase(k) = 2; % difficult+no unc / easy+uncertain
    elseif k > 30 && k <= 40
        phase(k) = 3; % difficult + uncertain
    elseif k > 50 && k <= 60
        phase(k) = 3; % difficult + uncertain
    elseif k > 60
        phase(k) = 2; % difficult+no unc / easy+uncertain
    end
end


%% create data struct
dat.trials = trials;
dat.phase = phase;
dat.u_bin = u_bin;
dat.y_pred = y_pred;
dat.y_mc = y_mc;
dat.y_c = y_c;
dat.y_tol = y_tol;
dat.y_av = y_av;
dat.u_pe = u_pe;
dat.u_pe_sq = u_pe_sq;
dat.u_pe_pos = u_pe_pos;
dat.u_pe_neg = u_pe_neg;
dat.u_mab2 = u_mab2;
dat.u_mab4 = u_mab4;
dat.exp_seq = exp_seq;


%% save struct
save('data\discovery_set_tmp.mat', 'dat', '-mat');


end