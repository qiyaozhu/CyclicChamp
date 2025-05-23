% This function generates cyclic peptide backbones with possible disulfide
% cross-links using simulated annealing.
% @ n: number of residues in the backbone
% @ N: number of simulated annealing trajectories to try
% @ best_k: for each trajectory, write out the best_k backbones found as pdbs
% @ outdir: the output directory
% @ is_disulfide: a boolean true/false to contain disulfide bond or not
% @ no_cys: if disulfide bond, positions to be excluded
% @ min_sep, max_sep: the min and max separation between the disulfide bond
% residues, bond (i, i+sep)
function Backbone_SA(n, N, best_k, outdir, is_disulfide, no_cys, min_sep, max_sep)

% Please add the correct paths to these two folders
addpath("/Users/qzhu/Desktop/CyclicChamp/my_functions");
addpath("/Users/qzhu/Desktop/CyclicChamp/RamaMap");

% Shuffle the random number generator
rng('shuffle');

% Ramachandran plot for symmetric glycine
Boundaries = load('ramabin_glycine.mat').Boundaries.';

% atom properties
[LJ_radius, LJ_well, ~, ~, ~, ~] = FA_parameter(n); % Lennard-Jones parameters
D = n_bonds(n); % number of bonds between any atom pair

% Read bond angle parameters for forming disulfide bond
[min_CaCbCb, max_CaCbCb, step_CaCbCb, data_CaCbCb] = read_histogram("centroid_CaCbCb_angle_score");
CaCbCb_x = [0, min_CaCbCb:step_CaCbCb:max_CaCbCb];
CaCbCb_y = [100; data_CaCbCb];

[min_CaCbCbCa, max_CaCbCbCa, step_CaCbCbCa, CaCbCbCa_y] = read_histogram("centroid_CaCbCbCa_dihedral_score");
CaCbCbCa_x = min_CaCbCbCa:step_CaCbCbCa:max_CaCbCbCa;

[min_NCaCaC, max_NCaCaC, step_NCaCaC, NCaCaC_y] = read_histogram("centroid_backbone_dihedral_score");
NCaCaC_x = min_NCaCaC:step_NCaCaC:max_NCaCaC;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% This is the section that you could play with to change 
%%%%%%%%%%%%%%%% parameters used for simulated annealing. I found often an
%%%%%%%%%%%%%%%% acceptance rate ~15% of random moves would work nicely.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% weights for different energy terms
w_cyc = 15;
w_rep = 0.4;
w_hbd_ramping = 10;
if is_disulfide
    w_disul = 1;
else
    w_disul = 0;
end

% parameters for simulated annealing
t0_ramping = w_cyc*5 + w_rep*n*10 + w_hbd_ramping*2 + w_disul*6; % initial temperature, the higher temperature, the more likely to accept bad random moves
c = 40; % temperature dropping rate
k0 = 0.7; % initial random move disc radius, the larger the radius, the larger the perturbations to the torsion angles
b = 16; % disc radius shrinking rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Parameters in this section are set to their default
%%%%%%%%%%%%%%%% values. No need to change unless you have different
%%%%%%%%%%%%%%%% requirements for the backbones.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 20000; % total number of steps in a simulated annealing trajectory
ramp0 = 1; % initial ramping factor, default is always 1, the most lenient requirements on H-bonds and disulfide bond

penalty = 2; % penalty score assigned to oversaturated hbonds to carbonyls
reward = 2; % reward score assigned to long-range hbonds or no floppy unbonded regions

% criteria for good backbone candidates
rep_cri = 10 + (n-7)*10/17;
cyc_cri = 0.1;
count_cri = ceil(n/3);
disul_cri = -6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Six torsion bin centers of the symmetric glycine Rama space
InitAngles = [  -2, -1.7, -1.9,  1.9,   1.7,   2; ...
              -2.4, 0.17,  2.6, -2.6, -0.17, 2.4];

% Randomly generate N initial backbones with torsion angles selected from
% the six torsion bin centers
IP = reshape(randsample(size(InitAngles, 2), n*N, true), [n, N]);
Angles_Initial = zeros(2*n, N);
for bb = 1 : N
    bin = IP(:, bb);
    Angles_Initial(:, bb) = reshape(InitAngles(:,bin),[],1);
end

% Record all the data for the good backbones found within each trajectory
% that satisfy the criteria
BEST_CAND_DATA = cell(1, N);
parfor repeat = 1 : N
    best_cand = [];
    best_score = [];
    best_data = [];
    best_coordinates = [];
    best_cys_res1 = [];
    best_cys_res2 = [];
    best_cys_name1 = [];
    best_cys_name2 = [];

    % initial coordinates and energies
    angles_initial = Angles_Initial(:,repeat);
    [coordinates, cyc] = peptide_cyc(angles_initial);
    rep = fa_rep_manhattan(coordinates, D, LJ_radius, LJ_well);
    [hbond_ramping, ~, ~] = E_hbond_ramping(coordinates, ramp0, penalty, reward);
    if is_disulfide
        [E_disul_ramping, ~, ~, ~, ~, ~] = backbone_disulfide(coordinates, n, no_cys, CaCbCb_x, CaCbCb_y, CaCbCbCa_x, CaCbCbCa_y, NCaCaC_x, NCaCaC_y, ramp0, min_sep, max_sep);
        E_total_ramping = w_cyc*cyc + w_rep*rep + w_hbd_ramping*hbond_ramping + w_disul*E_disul_ramping;
    else
        E_total_ramping = w_cyc*cyc + w_rep*rep + w_hbd_ramping*hbond_ramping;
    end

    % This records the number of random moves accepted, could be used to
    % adjust the parameters
    accept_ramping = 0;

    % Begin the simulated annealing trajectory
    angles = angles_initial;
    cand = 0;
    for i = 1 : M
        % generate new random move
        k = k0/(1+b*i/M);
        p = random_move(k,n);
        [~, angles_new] = feasible(p, angles, Boundaries);

        % gradually decrease the rampling/relaxation on the energy
        % calculation to go back to the original energy models
        ramp = ramp0*(M-i)/M;

        % compute energies for the new backbone after the random move
        [coordinates, cyc] = peptide_cyc(angles_new);
        rep = fa_rep_manhattan(coordinates, D, LJ_radius, LJ_well);
        [hbond_ramping, count, oversat] = E_hbond_ramping(coordinates, ramp, penalty, reward);
        if is_disulfide
            [E_disul_ramping, E_disul, cys_res1, cys_res2, cys_name1, cys_name2] = backbone_disulfide(coordinates, n, no_cys, CaCbCb_x, CaCbCb_y, CaCbCbCa_x, CaCbCbCa_y, NCaCaC_x, NCaCaC_y, ramp, min_sep, max_sep);
            E_total_ramping_new = w_cyc*cyc + w_rep*rep + w_hbd_ramping*hbond_ramping + w_disul*E_disul_ramping;
        else
            E_total_ramping_new = w_cyc*cyc + w_rep*rep + w_hbd_ramping*hbond_ramping;
        end
        
        % decide to accept the move or not using Metropolis criterion
        pass = false;
        temp = t0_ramping/(1+c*i/M);
        if E_total_ramping_new <= E_total_ramping
            pass = true;
        else
            prob = exp(1)^((E_total_ramping-E_total_ramping_new)/temp);
            if rand <= prob
                pass = true;
            end
        end

        % update the current status if move accepted
        if pass
            % display("i="+i+", total="+E_total_ramping_new+", cyc="+cyc+", rep="+rep+", hbond="+hbond_ramping+", hcount="+count+", disulfide="+E_disul_ramping);
            angles = angles_new;
            E_total_ramping = E_total_ramping_new;
            accept_ramping = accept_ramping + 1;

            % record the data if satisfy good backbone criteria
            satisfy = false;
            if is_disulfide
                if rep <= rep_cri && cyc <= cyc_cri && count >= count_cri && oversat <= 0 && E_disul <= disul_cri
                    satisfy = true;
                end
            else
                if rep <= rep_cri && cyc <= cyc_cri && count >= count_cri && oversat <= 0
                    satisfy = true;
                end
            end
            
            if satisfy
                cand = cand + 1;
                display("SA"+repeat+"_cand"+cand+", total="+E_total_ramping+", cyc="+cyc+", rep="+rep+", hcount="+count);
                best_cand = [best_cand, cand];
                best_score = [best_score, E_total_ramping];
                if is_disulfide
                    best_data = [best_data; cyc, rep, count, E_disul];
                else
                    best_data = [best_data; cyc, rep, count];
                end
                best_coordinates = [best_coordinates, coordinates];
                if is_disulfide
                    best_cys_res1 = [best_cys_res1, cys_res1];
                    best_cys_res2 = [best_cys_res2, cys_res2];
                    best_cys_name1 = [best_cys_name1, cys_name1];
                    best_cys_name2 = [best_cys_name2, cys_name2];
                end
            end
        end
    end

    % Output the good backbones as pdb files
    if ~isempty(best_score)
        [~, best_ind] = mink(best_score, best_k);
        best_data_rows = [];
        file_list = strings(length(best_ind), 1);

        for best_i = 1:length(best_ind)
            cand_ind = best_ind(best_i);
            filename = outdir+"SA"+repeat+"_cand"+best_cand(cand_ind)+".pdb";
            file_list(best_i) = filename;
            best_data_rows = [best_data_rows; best_data(cand_ind,:)];
            if is_disulfide
                plot_peptide_disulf(best_coordinates(:, cand_ind*3-2:cand_ind*3), n, filename, best_cys_res1(cand_ind), best_cys_res2(cand_ind), best_cys_name1(cand_ind), best_cys_name2(cand_ind));
            else
                plot_peptide(best_coordinates(:, cand_ind*3-2:cand_ind*3), n, filename);
            end
        end

        BEST_CAND_DATA{repeat}.filenames = file_list;
        BEST_CAND_DATA{repeat}.data = best_data_rows;
    end

    display("SA = "+repeat+", Ramp accept = "+accept_ramping);
end

save(outdir+"Backbone_data.mat", 'BEST_CAND_DATA');
end




% Read the histograms for disulfide energy calculation
function [min_val, max_val, step_size, data_values] = read_histogram(hist)
fid = fopen(hist, 'r');
min_val = NaN;
max_val = NaN;
step_size = NaN;
data_values = [];

while ~feof(fid)
    line = fgetl(fid);
    
    if contains(line, '@minimum')
        min_val = sscanf(line, '@minimum %f');
    elseif contains(line, '@maximum')
        max_val = sscanf(line, '@maximum %f');
    elseif contains(line, '@step')
        step_size = sscanf(line, '@step %f');
    elseif ~startsWith(line, '#') && ~startsWith(line, '@') && ~isempty(line)
        % If the line is not a comment and not empty, it’s a data line
        value = sscanf(line, '%f');
        data_values = [data_values; value];
    end
end

if length(data_values) ~= (max_val-min_val)/step_size
    error("histogram bin number is not correct!");
end
data_values = [data_values; data_values(end)];

fclose(fid);
end


% Plot the peptide backbone with disulfide bond
function plot_peptide_disulf(coordinates, n, filename, cys_res1, cys_res2, cys_name1, cys_name2)
peptide.X = coordinates(:,1).';
peptide.Y = coordinates(:,2).';
peptide.Z = coordinates(:,3).';

peptide.atomNum = 1 : length(peptide.X);
peptide.atomName = repmat(["N", "H", "CA", "1HA", "2HA", "C", "O"], 1, n);

pepres_names = repmat("GLY", 1, 7*n);
% Change the cysteine residue names
pepres_names(7*cys_res1-6 : 7*cys_res1) = repmat(cys_name1, 1, 7);
pepres_names(7*cys_res2-6 : 7*cys_res2) = repmat(cys_name2, 1, 7);
peptide.resName = pepres_names;

peptide.resNum = repelem(1:n, 7);
peptide.element = repmat(["N", "H", "C", "H", "H", "C", "O"], 1, n);

peptide.outfile = filename;
file = fopen(filename, "w");
fclose(file);
mat2pdb(peptide);
end


% Plot the peptide backbone
function plot_peptide(coordinates, n, filename)
peptide.X = coordinates(:,1).';
peptide.Y = coordinates(:,2).';
peptide.Z = coordinates(:,3).';

peptide.atomNum = 1 : length(peptide.X);
peptide.atomName = repmat(["N", "H", "CA", "1HA", "2HA", "C", "O"], 1, n);

peptide.resName = repmat("GLY", 1, 7*n);
peptide.resNum = repelem(1:n, 7);
peptide.element = repmat(["N", "H", "C", "H", "H", "C", "O"], 1, n);

peptide.outfile = filename;
file = fopen(filename, "w");
fclose(file);
mat2pdb(peptide);
end