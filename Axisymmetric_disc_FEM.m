%%----------------------------------------------------------------%%
%%%% Calculates displacement, stresses and strains for a
%%%% 2-dimensional spinning disc using Finite element method.
%%----------------------------------------------------------------%%

clear;
close all;

%% Choose constants

rho = 1; % density
omega = 1; % angular velocity
nu = 0.3; % poissons number
E_modul = 1; % E module

Ri = 1; % Inner radius
Ry = 6; % Outer radius

N_element = 16; % Number of elements

%% FEM Calculations

Constants = [Ri,Ry,N_element,rho,omega,nu,E_modul];

[K, F, nodes, dr] = Axisym_stiffness_load(Constants); % calculate stiffness and load vector

% Check condition number
cond_no=rcond(K);
cond_no=1./cond_no;
disp(['Condition number is ' num2str(cond_no)])

U = K\F; % displacement vektor

% Returns radial and theta stresses with elements and a radial vector for 
% the analytic solution and element positions
[r_analytic, sigmaR, sigmatheta, Sigma_r_element, Sigma_theta_element, rel_element] = Axisym_stresses(Constants,U,dr,nodes);


%% Creates functions for analytic solution

% Analytic solution - for stresses

A = (3+nu)/8 * rho * omega^2 * (Ry^2 + Ri^2);
B = (3+nu)/8 * rho * omega^2 * Ry^2 * Ri^2;
Cr = (3+nu)/8 * rho * omega^2;
Ctheta = (1+3*nu)/8 * rho * omega^2;

C = [A B Cr Ctheta];

% Analytic solution - for displacement vector

C3 = ((1-nu^2)/(8*E_modul)) * rho * omega^2;
C2 = (1+nu) * (3+nu) * rho * omega^2 * Ri^2 * Ry^2 /(8*E_modul);
C1 = (3+nu) * (1-nu) * rho * omega^2 * (Ry^2 + Ri^2) / (8*E_modul);

C0 = [C1 C2 C3];

if Ri == 0

    % Analytic displacement
    Analytic_U_func = @(r) C0(1).*r  - C0(3) * r.^3;
    % Radial stress 
    Analytic_sigmar_func = @(r) C(1) - C(3).*r.^2;
    % Theta stress
    Analytic_sigmatheta_func = @(r) C(1) - C(4).*r.^2;
    % epsilon strains
    epsilonr = @(r) C0(1) - 3 * C0(3) .* r.^2; % du/dr
    epsilontheta = @(r) C0(1) - C0(3) .* r.^2; % u/r

else

    % Analytic displacement
    Analytic_U_func = @(r) C0(1).*r + C0(2)./(r) - C0(3) * r.^3;
    % Radial stress 
    Analytic_sigmar_func = @(r) C(1) - C(2)./(r.^2) - C(3).*r.^2;
    % Theta stress
    Analytic_sigmatheta_func = @(r) C(1) + C(2)./(r.^2) - C(4).*r.^2;
    % epsilon strains
    epsilonr = @(r) C0(1) - C0(2)./(r.^2) - 3 * C0(3) .* r.^2; % du/dr
    epsilontheta = @(r) C0(1) + C0(2)./(r.^2) - C0(3) .* r.^2; % u/r

end    

% Creates vector for analytic displacement
Analytic_U = Analytic_U_func(r_analytic);

 % Sigmar plot
Analytic_sigmar = Analytic_sigmar_func(r_analytic); % Analytic sigma r
Analytic_sigmar_nodes = Analytic_sigmar_func(nodes); % for error control

Analytic_sigmatheta = Analytic_sigmatheta_func(r_analytic); % Analytic Sigma theta
Analytic_sigmatheta_nodes = Analytic_sigmatheta_func(nodes); % for error control

SigmaR_exact = Analytic_sigmar_func(rel_element); % Element wise sigmar analytic
Sigmatheta_exact = Analytic_sigmatheta_func(rel_element); % Element wise sigmatheta analytic

%% Calculates error for displacement and stresses

Analytic_U_nodes = Analytic_U_func(nodes);
Analytic_U_el_mid = Analytic_U_func(r_analytic);

U_interp = interp1(nodes,U',r_analytic); % interpolate 100 values in domain

% L_2 norm error for u(r)
L2_error_U_nodes = norm(Analytic_U_nodes - U');
L2_error_U_interp = norm(Analytic_U_el_mid - U_interp);
disp(['L2 norm error U with interpolation ' num2str(L2_error_U_interp) '']) % with interpolation
disp(['L2 norm error U at nodes ' num2str(L2_error_U_nodes) '']) % at nodes

% L_2 norm error for sigma-r 
L_2_sigmaR_elementwise = norm(SigmaR_exact-Sigma_r_element);
L_2_sigmaR = norm(Analytic_sigmar_nodes - sigmaR);
disp(['L2 norm error of sigma_r smoothed at nodes ' num2str(L_2_sigmaR) '']) % at nodes
disp(['L2 norm error of sigma_r elemenwise at nodes ' num2str(L_2_sigmaR_elementwise) '']) % at nodes

% L_2 norm error for sigma-theta
L_2_sigmatheta_elementwise = norm(Sigmatheta_exact-Sigma_theta_element);
L_2_sigmatheta = norm(Analytic_sigmatheta_nodes - sigmatheta);
disp(['L2 norm error of sigma_theta smoothed at nodes ' num2str(L_2_sigmatheta) '']) % at nodes
disp(['L2 norm error of sigma_theta elementwise at nodes ' num2str(L_2_sigmatheta_elementwise) '']) % at nodes

%% Finds energynorm and calculates error

% Energy norm FEM
energy_norm_FEM = 0.5 * U' * K * U - U' * F; 
disp(['Energynorm FEM - case ' num2str(energy_norm_FEM)])

% Energy norm Analytic
integral1 = @(r) 0.5 .* epsilonr(r) .* Analytic_sigmar_func(r) .* r + 0.5 .* epsilontheta(r) .* Analytic_sigmatheta_func(r) .* r;
Q1 = integral(@(x) integral1(x), Ri, Ry);

integral2 = @(r) Analytic_U_func(r) .* rho * omega^2 .* r.^2;
Q2 = integral (@(x) integral2(x), Ri, Ry);

energy_norm_analytic = Q1 - Q2;
disp(['Energynorm Analytic - case ' num2str(energy_norm_analytic)])

%Energynorm difference FEM - Analytic
energy_norm_skillnad = energy_norm_FEM - energy_norm_analytic;
disp(['Energynorm difference FEM - Analytic ' num2str(energy_norm_skillnad)])

%% Creates plots of results

% All FEM calculations in one graph
figure(1)
plot(nodes, U);
hold on
plot(nodes,sigmatheta)
plot(nodes,sigmaR)
ylabel('Amplitude');
xlabel('Radial distance r')
legend('U(r)','sigmatheta(r)','sigmar(r)')
hold off

figure(2)
plot(r_analytic,Analytic_U);
hold on
plot(nodes,U,'-x');
xlabel('Radial distance r');
ylabel('Displacement U(r)')
legend('Analytic-sol', 'FEM-sol')
hold off

figure(3)
plot(r_analytic,Analytic_sigmar);
hold on
plot(nodes,sigmaR,'-x');
plot(rel_element,Sigma_r_element,'-x')
ylabel('Radial stress Sigma-r');
xlabel('Radial distance r');
legend('Analytic-sol','FEM-sol','Element');
hold off

figure(4)
plot(r_analytic,Analytic_sigmatheta);
hold on
plot(nodes,sigmatheta,'-x');
plot(rel_element,Sigma_theta_element,'-x');
ylabel('Radial stress Sigma-theta');
xlabel('Radial distance r');
legend('Analytic-sol','FEM-sol','Element');
hold off
