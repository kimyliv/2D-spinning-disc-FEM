function [K, F, nodes, dr] = Axisym_stiffness_load(Constants)


Ri = Constants(1);
Ry = Constants(2);
N_element = Constants(3);
rho = Constants(4);
omega = Constants(5);
nu = Constants(6);
E_modul = Constants(7);


dr = (Ry - Ri) / N_element; % Length each element
nodes = linspace(Ri, Ry, N_element + 1); % Node positions

% Gauss points och weights
gauss_points = [-1/sqrt(3) 1/sqrt(3)];
% gauss_punkter = [-sqrt(3/5) 0 sqrt(3/5)];
% gauss_punkter = [0];
weights = [1 1];
% vikter = [5/9 8/9 5/9];
% vikter = [2];
% gauss_punkter = [-sqrt(5 + 2*sqrt(10/7))/3, -sqrt(5 - 2*sqrt(10/7))/3, 0, sqrt(5 - 2*sqrt(10/7))/3, sqrt(5 + 2*sqrt(10/7))/3];
% vikter = [(322 - 13*sqrt(70))/900, (322 + 13*sqrt(70))/900, 128/225, (322 + 13*sqrt(70))/900, (322 - 13*sqrt(70))/900];


K = zeros(N_element + 1); % Global Stiffness
F = zeros(N_element + 1, 1); % Global Force

for k = 1:N_element

    k_local = zeros(2,2); % Lokal Stiffness
    f_local = zeros(2,1); % Lokal Force

    for i = 1:length(weights)

        s = gauss_points(i);
        w = weights(i);

        N = [(1-s)/2 (1+s)/2]; % 2-Node functions
        
        r_mid = (nodes(k+1) + nodes(k))/2;
        r = r_mid + s*dr/2;

        J = dr/2; % Jacobian

        B = [-1/dr 1/dr;N(1)/r  N(2)/r];
        E = [1 nu; nu 1];
        
        k_local = k_local + w * J * (E_modul/(1-nu^2)) * (B' * (E * (B * r)));
        f_local = f_local + w * rho * omega^2 * J * r^2 * N';

    end

K(k:k+1,k:k+1) = K(k:k+1,k:k+1) + k_local;
F(k:k+1,1) = F(k:k+1,1) + f_local;

end
