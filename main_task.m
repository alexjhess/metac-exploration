%% [METAC] exploration data set - task data


%% remove default toolboxes that could be conflicting

% List all paths in search path and break them into a cell array
dirs = regexp(path,['[^;]*'],'match');

% Return all search path directories containing this string
strsToFind = {'tapas'}; 

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


%% load data
% dat = load_discovery_set();
load('data\discovery_set_tmp.mat');


%% visualise raw data
metac_plot_raw_data(dat);




%% MC: cont eHGF on raw PE (scale of PE ???)
dat = metac_raw_PE_mc_ehgf(dat);

%% INT: 3-level binary HGF MAB
dat = metac_int_hgf_binary_mab(dat);

%% INT: 3-level binary eHGF MAB (ERROR..)
dat = metac_int_ehgf_binary_mab(dat);


%% INT: classical 3-level binary eHGF (NOT APPROPRIATE!!!)
dat = metac_int_ehgf_binary(dat);






%% ____________
% build autoregr mc obs model (giuliara). HGF MAB ???

dat.ehgf_mc_pe.sub(n).sim = tapas_simModel(dat.u_bin(:,n),...
    'tapas_ehgf_binary',...
    [NaN 0 1 NaN 1 1 NaN 0 0 1 1 NaN -2 2],...
    'mc_autoreg_pe_inobs',... 'mc_null_inobs',...
    [0.05 1 -2 -1.5 -0.04 0.005],... [0.05 -2 0.005],...
    12345 ...
    );
tapas_ehgf_binary_plotTraj(dat.ehgf_mc_pe.sub(n).sim)
hold on
plot(dat.ehgf_mc_pe.sub(n).sim.y(:,2), '.')

%% sim null model. HGF MAB ???
dat.ehgf_mc_null.sub(n).sim = tapas_simModel(dat.u_bin(:,n),...
    'tapas_ehgf_binary',...
    [NaN 0 1 NaN 1 1 NaN 0 0 1 1 NaN -2 2],...
    'mc_null_inobs',...
    [0.05 0.5 0.005],...
    12345 ...
    );
tapas_ehgf_binary_plotTraj(dat.ehgf_mc_null.sub(n).sim)
hold on
plot(dat.ehgf_mc_null.sub(n).sim.y(:,2), '.')


%% fit metacognitive HGF (giuliara). HGF MAB ???

for n = 1:size(dat.u_bin,2)
    dat.ehgf_mc_null.sub(n).est = tapas_fitModel([dat.y_pred(:,n), dat.y_mc(:,n)],...
        dat.u_bin(:,n),...
        tapas_ehgf_binary_config,...
        mc_null_inobs_config,...
        tapas_quasinewton_optim_config ...
        );
    tapas_ehgf_binary_plotTraj(dat.ehgf_mc_null.sub(n).est)
    hold on
    plot(dat.ehgf_mc_null.sub(n).est.y(:,2), '.')
end




