% Function to generate random move direction
% @ k: the move is within a disc of radius k and uniformly distributed
% @ n: the number of residues
% @ p: a vector of length 2n, contains the perturbation/move to be added to
% phi_1, psi_1, phi_2, psi_2, ..., phi_n, psi_n
function p = random_move(k, n)

% the angle direction within the disc
a = 2*pi*rand(1,n);

% the step size towards the angle direction
% note that the square root is needed to achieve uniform distribution, 
% otherwise the points will accumulate near the center
r = k*sqrt(rand(1,n));

p = reshape([r.*sin(a); r.*cos(a)], [2*n,1]);

end