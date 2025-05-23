% Function to find a feasible direction within the Rama space
% @ p: initially proposed random moves for n residues
% @ angles: current torsion angles phi_1, psi_1, ..., phi_n, psi_n
% @ Boundaries: Ramachandran space boundary vertices, polygon assumed
% @ p_f: feasible move within the Rama space
% @ angles_f: angles after the feasible move
function [p_f, angles_f] = feasible(p, angles, Boundaries)
p_f = p;
angles_f = angles + p_f;
angles_f = mod(angles_f+pi, 2*pi) - pi;

% Check after the move, which residues' phi, psi are still within the Rama space
in_or_not = double(inpolygon(angles_f(1:2:end-1), angles_f(2:2:end), Boundaries(1,:), Boundaries(2,:)).');
in_or_not = repelem(in_or_not, 2);

% Try retract moves by a factor of 0.8 in the opposite direction if go outside of the Rama space
while sum(in_or_not==0)>0
    in_or_not(in_or_not==0) = -0.8;
    p_f = p_f.*in_or_not.';
    angles_f = angles + p_f;
    angles_f = mod(angles_f+pi, 2*pi) - pi;
    in_or_not = double(inpolygon(angles_f(1:2:end-1), angles_f(2:2:end), Boundaries(1,:), Boundaries(2,:)).');
    in_or_not = repelem(in_or_not, 2);
end
end