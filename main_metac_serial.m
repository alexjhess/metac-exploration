%
% Main script for running the analysis of the METAC 1 study locally.
%
%
% _________________________________________________________________________
% Author: Alex Hess
%
% Copyright (C) 2025 Translational Neuromodeling Unit
%                    Institute for Biomedical Engineering
%                    University of Zurich & ETH Zurich
%
% This file is released under the terms of the GNU General Public Licence
% (GPL), version 3. You can redistribute it and/or modify it under the
% terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any
% later version.
%
% This file is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
% more details.
% 
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <https://www.gnu.org/licenses/>.
% _________________________________________________________________________

%% specify the number of local cores (for MATLAB's parallel toolbox)
local_cores = 4;

%% set global variables
EULER = 0;
nMod = 5;
nS_ds = 20;
nS_vs = 30;
nS_sim = 100;

%% remove default toolboxes that could be conflicting

% List all paths in search path and break them into a cell array
dirs = regexp(path,['[^;]*'],'match');

% Return all search path directories containing this string
strsToFind = {'tapas', 'hgf-toolbox', 'comb_obs_models', 'mc_obs', 'VBA-toolbox'}; 

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
addpath('mc_obs');
cd('VBA-toolbox')
VBA_setup();
cd ..

%% ________________________________________________________________________
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DISCOVERY SET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% step 1:
% set specs, load pilot data, model-agnostic analysis, specify model space
metac_job_1_init(EULER);

%% step 2:
% fit behavioural models on discovery set
% request cores
local = parcluster('local'); local.NumWorkers = local_cores;
pool = parpool(local, local.NumWorkers);

% model inversion
parfor n = 1:nS_ds
    for m = 1:nMod
        metac_job_2_ds_fit(EULER, n, m)
    end 
end

% shuting down cores
pool.delete()

disp('model inversion on the discovery set successfully completed')

%% step 3:
metac_job_3_ep(EULER);

%% ________________________________________________________________________
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% step 4:
% setup simulation analyses (individual model level)
metac_job_4_sim(EULER);

%% step 5:
% inversion of synthetic data sets

% request cores
local = parcluster('local');
local.NumWorkers = local_cores;
pool = parpool(local, local.NumWorkers);

% model inversion
parfor n = 1:nS_sim
    for m = 1:nMod
        for i = 1:nMod
            metac_job_5_sim_data_fit(EULER, n, m, i)
        end
    end 
end

% shuting down cores
pool.delete()

disp('model inversion on the synthetic data set successfully completed')

% %% step 6:
% % recovery analysis (individual model level)
% job_runner_6_rec_analysis(EULER);
% 
% %% step 7:
% % family-level simulations
% job_runner_7_family_power(EULER);
% 
% %% step 8:
% % perform family-level BMS on synthetic data sets
% 
% % request cores
% local = parcluster('local'); local.NumWorkers = local_cores;
% pool = parpool(local, local.NumWorkers);
% 
% % run family level BMS for all different fam sim opts
% parfor j = 1:nSim_fam
%     job_runner_8_sim_family_BMS(EULER, j)
% end
% 
% disp('family-level simulations completed')
% 
% % shuting down cores
% pool.delete()
% 
% %% step 9:
% % analyze results from family BMS sim
% job_runner_9_family_power_res(EULER)
% 
% %% ________________________________________________________________________
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN DATA SET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %% step 10:
% % load main data set, model-agnostic analysis
% job_runner_10_main_model_agnostic(EULER);
% 
% %% step 11:
% % fit main data set
% 
% % model inversion
% parfor n = 1:nS_vs
%     for m = 1:nMod
%         job_runner_11_main_modinv(EULER, n, m)
%     end 
% end
% 
% % shuting down cores
% pool.delete()
% 
% disp('model inversion on the main data set successfully completed')
% 
% %% step 12:
% % model-based analysis of main data set
% job_runner_12_main_results(EULER);
% 
% %% ________________________________________________________________________
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%% OPTIONAL PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %% [optional] step 13:
% % create figures
% job_runner_plot_results(EULER);
% 
% %% [optional] step 14:
% % print stats that are reported in the paper
% job_runner_print_stats(EULER);
% 
% %% [optional] step 15:
% % create figures for paper
% job_runner_paper_figs(EULER)
% 
% %% [optional] step 16:
% % create figures for CCN Conference paper (2023)
% job_runner_ccn_figs(EULER);
