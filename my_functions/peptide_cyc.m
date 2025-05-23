%% Glycine each residue has 7 atoms: N, H, CA, 1HA, 2HA, C, O
% @ angles: an input vector containing phi_1, psi_1, phi_2, psi_2, ..., phi_n, psi_n
% @ coordinates: the output atom coordinate matrix of size 7n x 3
% @ error: the cyclic error measuring how far this backbone is away from a perfect closure
% Note: To compute the cyclic error, we add a virtual atom C before the N
% terminus, and two virtual atoms N and CA after the C terminus.

function [coordinates, error] = peptide_cyc(angles)
% bond angles
COR = deg2rad(120.8);
NHR = deg2rad(119.2);

% bond lengths
NCa = 1.458;
CaH = 1.090;
CaC = 1.524;
CN = 1.329;
CO = 1.231;
NH = 1.010;

% pre-computed rotation matrices caused by the bond angles and omega torsion
T_NR = [0.5255, -0.8508, 0; 0.8508, 0.5255, 0; 0, 0, 1];
T_CaR = [0.3616, -0.9323, 0; 0.9323, 0.3616, 0; 0, 0, 1];
CR_omega = [0.4415, 0.8973, 0; 0.8973, -0.4415, 0; 0, 0, -1];

% pre-computed translation vectors for adding the H and O atoms
NH_trans = [NH*cos(pi-NHR); -NH*sin(pi-NHR); 0];
CO_trans = [CO*cos(pi-COR); -CO*sin(pi-COR); 0];

% initialize variables
n = length(angles)/2;
coordinates = zeros(7*n, 3);
C = [0;0;0]; % start from virtual C at the origin
rotation = eye(3);

% get the backbone coordinates
for i = 1 : n+1
    if i <= n
        N = C + rotation*[CN; 0; 0];
        coordinates(7*i-6,:) = N.';

        H = N + rotation*NH_trans;
        coordinates(7*i-5,:) = H.';

        rotation = rotation*T_NR*R(angles(2*i-1));
        Ca = N + rotation*[NCa; 0; 0];
        coordinates(7*i-4,:) = Ca.';

        if i == n
            rotation_O = rotation;
        end

        rotation = rotation*T_CaR*R(angles(2*i));
        C = Ca + rotation*[CaC; 0; 0];
        coordinates(7*i-1,:) = C.';

        % The last O needs to be fixed by computing the actual psi_n
        if i < n
            O = C + rotation*CO_trans;
        else
            last_psi = torsion(N, Ca, C, coordinates(1,:).');
            O = C + rotation_O*T_CaR*R(last_psi)*CO_trans;
        end
        coordinates(7*i,:) = O.';

        rotation = rotation*CR_omega;

    else
        virtual_N = C + rotation*[CN; 0; 0];
        rotation = rotation*T_NR;
        virtual_Ca = virtual_N + rotation*[NCa; 0; 0];
    end
end

% To fix the first H, need to go in reverse order, C1-CA1-N1-C_last
x = coordinates(3,:)-coordinates(6,:);
x = x/norm(x);
u = coordinates(1,:) - coordinates(6,:);
y = u - dot(u,x)*x;
y = y/norm(y);
z = cross(x, y);

phi_first = torsion(coordinates(6,:), coordinates(3,:), coordinates(1,:), coordinates(7*n-1,:));
t = T_CaR*R(phi_first)*NH_trans;
H = coordinates(1,:) + t(1)*x + t(2)*y + t(3)*z;
coordinates(2,:) = H.';

% add the two hydrogen atoms at each Ca
for i = 1 : n
    Ncoord = coordinates(7*i-6,:);
    CAcoord = coordinates(7*i-4,:);
    Ccoord = coordinates(7*i-1,:);
    hor = -(Ncoord + (Ccoord-Ncoord)*1.20315/2.46077 - CAcoord);
    horcomp = hor / norm(hor) * CaH*cos(0.938693);
    pep = cross(Ccoord-CAcoord, Ncoord-CAcoord);
    pepcomp = pep / norm(pep) * CaH*sin(0.938693);
    coordinates(7*i-3,:) = CAcoord + horcomp + pepcomp;
    coordinates(7*i-2,:) = CAcoord + horcomp - pepcomp;
end

% Compute the cyclic error
x = virtual_N - coordinates(7*n-1,:).';
x = x/norm(x);
u = virtual_Ca - coordinates(7*n-1,:).';
y = u - dot(u,x)*x;
y = y/norm(y);
error = norm(coordinates(7*n-1,:))^2 + norm(x-[1;0;0])^2 + norm(y-[0;1;0])^2;
end


% Helper function for computing rotation matrix caused by torsion angles
function r = R(angle)
r = [1, 0, 0; 0, cos(angle), -sin(angle); 0, sin(angle), cos(angle)];
end

% Helper function for computing the torsion angle between four points
function chi = torsion(p1, p2, p3, p4)
b1 = p2-p1;
b2 = p3-p2;
b3 = p4-p3;
n1 = cross(b1, b2)/norm(cross(b1, b2));
n2 = cross(b2, b3)/norm(cross(b2, b3));
x = dot(n1, n2);
y = dot(cross(n1, n2), b2/norm(b2));
chi = atan2(y, x);
end