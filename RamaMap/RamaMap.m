%% Read Rosetta Glycine rama map
clc;
clear all;

file = fopen('ala.dat','rt');
C = textscan(file, '%s%f%f%f%f');
fclose(file);

Phi = C{2};
Psi = C{3};
Prob = C{4};
Energy = C{5};

F = reshape(Energy, 36, 36);
F = [F, F(:,1)];
F = [F; F(1,:)];
Energy = reshape(F, 37*37, 1);

P1 = reshape(Phi, 36, 36);
P1 = [P1, P1(:,1)];
P1 = [P1; -P1(1,:)];
Phi = reshape(P1, 37*37, 1);

P2 = reshape(Psi, 36, 36);
P2 = [P2, -P2(:,1)];
P2 = [P2; P2(1,:)];
Psi = reshape(P2, 37*37, 1);

% only need regions with energy<=10
goodIndices = find(Energy<=10);
phi = Phi(goodIndices);
psi = Psi(goodIndices);
energy = goodIndices;

save('ramabin_ala.mat', 'phi', 'psi', 'energy');

figure;
set(gcf,'color','w');
colormap autumn;
pointsize = 10;
scatter(phi, psi, pointsize, energy, 'filled');
colorbar;
xlim([-pi, pi]);
ylim([-pi, pi]);
xticks([-pi -2*pi/3 -pi/3 0 pi/3 2*pi/3 pi]);
yticks([-pi -2*pi/3 -pi/3 0 pi/3 2*pi/3 pi]);
xlabel('$\phi$','interpreter','latex');
ylabel('$\psi$','interpreter','latex');
set(gca,'FontSize',16,'FontWeight','Bold');
