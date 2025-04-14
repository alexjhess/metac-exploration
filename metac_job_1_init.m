function [] = metac_job_1_init(EULER)
% [] = metac_job_1_init(EULER)
%
% Init settings and paths for analysis.
%
% INPUT
%   EULER        binary           Binary indicator variable
%
%   OPTIONAL:
%
% OUTPUT    
%   argout       type
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

%% setup path
sdir = metac_paths(EULER);

%% analysis specifications
opts = metac_opts();

%% create results folder
metac_mkdir(opts,sdir);

%% load task data (EXPLORATION)
% dataset = 1; % 1=discovery set
% dat = load_task_data(dataset);
load(fullfile('data', opts.ds.fname));
ds.dat = dat;

%% preprocess task data
ds.pdat = metac_preproc(ds.dat);

%% create model space
[ds.ip.mod, bo] = metac_create_model_space(1); % 1=logit space


%% save model space and options as struct

% options
save_path1 = fullfile(sdir, 'results', ['options']);
save(save_path1, '-struct', 'opts');

% ds mod space
save_path = fullfile(sdir, 'results', 'ds', ['init_mod_space']);
save(save_path, '-struct', 'ds');

disp('model space ready for inversion on the discovery set.')

end
