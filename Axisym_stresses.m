function [r_analytic, sigmaR, sigmatheta,Sigma_r_element,Sigma_theta_element,rel_element] = Axisym_stresses(Constants,U,dr,nodes)

Ri = Constants(1);
Ry = Constants(2);
N_element = Constants(3);
nu = Constants(6);
E_modul = Constants(7);

sigmaR = zeros(1,N_element+1);
sigmatheta = zeros(1,N_element+1);
Sigma_r_element = zeros(1,2*N_element);
Sigma_theta_element = zeros(1,2*N_element);
rel_element = zeros(1,2*N_element);

rel_mid = linspace(Ri+dr/2,Ry-dr/2, N_element); % Radial value in middle of elements

for index = 1:N_element

    Umid = (U(index+1)+U(index))/2; % interpolates displacement between elements

    sigmavar_r = (U(index+1)-U(index))/dr + nu*(Umid/rel_mid(index));
    sigmaR(1,index) = sigmaR(1,index) + 0.5*sigmavar_r;
    sigmaR(1,index + 1) = sigmaR(1,index +1) + 0.5*sigmavar_r;
    Sigma_r_element(index*2-1) = sigmavar_r;
    Sigma_r_element(index*2) = sigmavar_r;

    sigmavar_theta = nu*(U(index+1)-U(index))/dr + (Umid/rel_mid(index));
    sigmatheta(1,index) = sigmatheta(1,index) + 0.5*sigmavar_theta;
    sigmatheta(1,index+1) = sigmatheta(1,index+1) + 0.5*sigmavar_theta;
    Sigma_theta_element(index*2-1) = sigmavar_theta;
    Sigma_theta_element(index*2) = sigmavar_theta;
    
    rel_element(2*index-1) = nodes(index);
    rel_element(2*index) = nodes(index+1);

end

sigmaR(1) = sigmaR(1)*2;
sigmaR(N_element + 1) = sigmaR(N_element + 1)*2;
sigmatheta(1) = sigmatheta(1)*2;
sigmatheta(N_element + 1) = sigmatheta(N_element + 1)*2;

sigmaR = E_modul/(1-nu^2).*sigmaR;
sigmatheta = E_modul/(1-nu^2).*sigmatheta;
Sigma_r_element = E_modul/(1-nu^2).*Sigma_r_element;
Sigma_theta_element = E_modul/(1-nu^2).*Sigma_theta_element;

r_analytic = linspace(Ri,Ry); % radial value for analytic solution
