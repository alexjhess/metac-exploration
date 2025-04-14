function [] = metac_mkdir(opts,sdir)
% [] = metac_mkdir(opts,sdir)
%
% Helper function to create directories for saving results & figures.
%
% INPUT
%   opts         struct        Settings for analysis pipeline
%   sdir         character (optional) Directory, where data should be stored.
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

if nargin < 2 || isempty(sdir)
    sdir =  pwd;
end

% create dir to store results
for n = 1:opts.ds.nS
    mkdir(fullfile(sdir, 'results', 'ds', ['sub', num2str(n)]));
end
for n = 1:opts.sim.n_sim
    mkdir(fullfile(sdir, 'results', 'sim', ['sub', num2str(n)]));
end

for n = 1:opts.vs.nS
    mkdir(fullfile(sdir, 'results', 'vs', ['sub', num2str(n)]));
end

mkdir(fullfile(sdir, 'results', 'sq2'));

% create dir to store figures
mkdir(fullfile(sdir, 'figures', 'ds', 'ip'))

mkdir(fullfile(sdir, 'figures', 'sim', 'traj'))
mkdir(fullfile(sdir, 'figures', 'sim', 'rec'))

mkdir(fullfile(sdir, 'figures', 'vs'))

end
