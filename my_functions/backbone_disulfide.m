% Reimplementation of Rosetta's centroid disulfide potential model with
% modifications. The original file is
% /main/source/src/core/scoring/disulfides/CentroidDisulfidePotential.cc
%
% A function to compute disulfide energy for all possible residue pairs in
% a backbone, trying both L- and D-CYS.
%
% @ backbone: coordinates of the backbone atoms (N, H, CA, 1HA, 2HA, C, O)
% @ n: size of peptide
% @ no_cys: the residues to be excluded from forming disulfide bond
% @ CaCbCb_x, CaCbCb_y: the data read from the histogram for Ca-Cb-Cb bond
% angle scores, use these for piecewise cubic Hermite interpolation
% @ CaCbCbCa_x, CaCbCbCa_y: the data read from the histogram for
% Ca-Cb-Cb-Ca dihedral score for piecewise cubic Hermite interpolation
% @ NCaCaC_x, NCaCaC_y: the data read from the histogram for N-Ca-Ca-C
% dihedral score for piecewise cubic Hermite interpolation
% @ ramp: ramp=1 most relaxed requirement for disulfide bond, ramp=0 normal
% requirement
% @ min_sep, max_sep: the min and max separation between the disulfide bond
% residues, bond (i, i+sep)
%
% @ score_ramp, score: the lowest disulfide energy from all possible pairs
% @ cys_res1, cys_res2: the lowest-energy pair residue numbers
% @ cys_name1, cys_name2: the lowest-energy pair residue names, L or D
function [score_ramp, score, cys_res1, cys_res2, cys_name1, cys_name2] = backbone_disulfide(backbone, n, no_cys, CaCbCb_x, CaCbCb_y, CaCbCbCa_x, CaCbCbCa_y, NCaCaC_x, NCaCaC_y, ramp, min_sep, max_sep)

score_ramp = 1000;
score = 1000;

for sep = min_sep : max_sep
    cys1 = 1 : n;
    cys2 = cys1 + sep;
    cys2(cys2>n) = cys2(cys2>n) - n;
    pairs = [cys1; cys2];
    pairs = sort(pairs, 1);
    pairs = unique(pairs.', "rows");

    for p = 1 : size(pairs, 1)
        c1 = pairs(p,1);
        c2 = pairs(p,2);

        if ~ismember(c1, no_cys) && ~ismember(c2, no_cys)
            % Calculate CB atoms
            N1 = backbone(7*c1-6,:);
            CA1 = backbone(7*c1-4,:);
            C1 = backbone(7*c1-1,:);
            N2 = backbone(7*c2-6,:);
            CA2 = backbone(7*c2-4,:);
            C2 = backbone(7*c2-1,:);

            for cn1 = ["CYS", "DCS"]
                for cn2 = ["CYS", "DCS"]
                    CB1 = compute_CB(N1, CA1, C1, lower(char(cn1)));
                    CB2 = compute_CB(N2, CA2, C2, lower(char(cn2)));

                    % Cb-Cb distance score, a weighted sum of 3 gaussians
                    Cb_Cb_dist = norm(CB1-CB2);
                    score_ramp_temp = cbcb_dist_score_ramping(Cb_Cb_dist^2, ramp);
                    score_temp = cbcb_dist_score(Cb_Cb_dist^2);

                    % Compute angle energy only if Cb-Cb distance is in range
                    if Cb_Cb_dist^2 > 10 && Cb_Cb_dist^2 < 400
                        % Use Cb-Cb distance to get factor for all other angle scores
                        score_factor = 0;
                        if Cb_Cb_dist^2 > 10 && Cb_Cb_dist^2 < 22
                            score_factor = 1;
                        elseif  Cb_Cb_dist^2 >= 22 && Cb_Cb_dist^2 < 400
                            score_factor = (Cb_Cb_dist^2-400) / (22-400);
                        end

                        % Get the CaCbCb bond angle component
                        CaCbCb_angle1 = acos(dot(CA1-CB1, CB2-CB1) / norm(CA1-CB1) / norm(CB2-CB1));
                        CaCbCb_angle2 = acos(dot(CA2-CB2, CB1-CB2) / norm(CA2-CB2) / norm(CB1-CB2));

                        % Angle scores are interpolated from histograms
                        angle_score = mean(score_factor * pchip(CaCbCb_x, CaCbCb_y, rad2deg([CaCbCb_angle1, CaCbCb_angle2])));
                        score_ramp_temp = score_ramp_temp + angle_score;
                        score_temp = score_temp + angle_score;

                        % Proceed to the most expensive dihedral angle energy
                        % calculations only if the score has chance of being
                        % lower than the current minimum
                        if score_ramp_temp < score_ramp + score_factor*3
                            CaCbCbCa_dihedral = -torsion(CA1, CB1, CB2, CA2);
                            if CaCbCbCa_dihedral < 0
                                CaCbCbCa_dihedral = CaCbCbCa_dihedral + 2*pi;
                            end

                            % For the backbone dihedral, needs to check L- or D-conformation
                            res_name1 = char(cn1);
                            res_name2 = char(cn2);
                            if res_name1(1) == "D"
                                if res_name2(1) == "D"
                                    NCaCaC_dihedral = -torsion(C1, CA1, CA2, N2);
                                else
                                    NCaCaC_dihedral = -torsion(C1, CA1, CA2, C2);
                                end
                            else
                                if res_name2(1) == "D"
                                    NCaCaC_dihedral = -torsion(N1, CA1, CA2, N2);
                                else
                                    NCaCaC_dihedral = -torsion(N1, CA1, CA2, C2);
                                end
                            end
                            if NCaCaC_dihedral < 0
                                NCaCaC_dihedral = NCaCaC_dihedral + 2*pi;
                            end

                            % Angle scores are interpolated from histograms
                            angle_score = score_factor * pchip(CaCbCbCa_x, CaCbCbCa_y, rad2deg(CaCbCbCa_dihedral)) ...
                                + score_factor * pchip(NCaCaC_x, NCaCaC_y, rad2deg(NCaCaC_dihedral));

                            % Cb-Cb distance score, a weighted sum of 3 gaussians
                            score_ramp_temp = score_ramp_temp + angle_score;
                            score_temp = score_temp + angle_score;
                        end
                    end

                    if score_ramp_temp < score_ramp
                        score_ramp = score_ramp_temp;
                        score = score_temp;
                        cys_res1 = c1;
                        cys_res2 = c2;
                        cys_name1 = cn1;
                        cys_name2 = cn2;
                    end
                end
            end
        end
    end
end
end



% Weighted sum of 3 gaussians for Cb-Cb distance score
function score = cbcb_dist_score(d)
means = [12.445, 15.327, 14];
sds = [1.1737973, 2.1955666, 0.3535534];
weights = [10.8864116, 33.5711622, 0.2658681];
score = -sum(weights ./ sds/sqrt(2*pi) .* exp(-(d-means).^2./(2*sds.^2)));
end


% Weighted sum of 3 gaussians for Cb-Cb distance score, with ramping
function score = cbcb_dist_score_ramping(d, ramp)
means = [12.445, 15.327, 14];
sds = [1.1737973, 2.1955666, 0.3535534];
weights = [10.8864116, 33.5711622, 0.2658681];

upper = 100;

if d <= 13.8364
    x = 13.8364 - (13.8364-d) / (13.8-6.5*(1-ramp)) * (13.8364-6.5);
else
    x = 13.8364 + (d-13.8364) / (22+(upper-22)*ramp - 13.8364) * (22-13.8364);
end

score = -sum(weights ./ sds/sqrt(2*pi) .* exp(-(x-means).^2./(2*sds.^2)));
end


% Function to compute Cb atom positions for both L- and D-CYS
function CB = compute_CB(N, CA, C, cys_name)
x = N - C;
x = x/norm(x);
u = CA - C;
y = u - dot(u,x)*x;
y = y/norm(y);
z = cross(x, y);

CNCaR = deg2rad(35.2504);
NCaCbR = deg2rad(110.6);
CaCbD = 1.528861;
if cys_name(1) == "d"
    CNCaCbT = deg2rad(121.6);
else
    CNCaCbT = deg2rad(-121.6);
end

rotation = T(CNCaR) * R(CNCaCbT) * T(NCaCbR);

t = rotation*[CaCbD; 0; 0];
CB = CA + t(1)*x + t(2)*y + t(3)*z;
end


function t = T(angle)
t = [cos(pi-angle), -sin(pi-angle), 0; sin(pi-angle), cos(pi-angle), 0; 0, 0, 1];
end

function r = R(angle)
r = [1, 0, 0; 0, cos(angle), -sin(angle); 0, sin(angle), cos(angle)];
end

function chi = torsion(p1, p2, p3, p4)
b1 = p2-p1;
b2 = p3-p2;
b3 = p4-p3;
n1 = cross(b1, b2)/norm(cross(b1, b2));
n2 = cross(b2, b3)/norm(cross(b2, b3));
m = cross(n1, b2/norm(b2));
x = dot(n1, n2);
y = dot(m, n2);
chi = atan2(y, x);
end