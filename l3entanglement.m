% Solving For L=N=3 entanglement
% Spatial Entanglement----------------------------------------------------------------------------------------
%% Finding the ground state eigen vector of 3x3 Hamiltonian 
clc;clear all
syms J U W real
sqrt12 = sqrt(12);
sqrt6 = sqrt(6);
% Define the 3x3 Hamiltonian Matrix
H = [0, -sqrt12, 0;
    -sqrt12, -3 + U, -sqrt6;
    0, -sqrt6, 3*U + W];

% Compute eigenvalues and eigenvectors
[V, D] = eig(H);

% Extract eigenvalues
eigenvalues = diag(D);

% Display eigenvalues
disp('Eigenvalues:');
disp(eigenvalues);

% Display eigenvectors
disp('Eigenvectors:');
disp(V);

% Display the third eigenvector
disp('Third Eigenvector:');
disp(V(:,3)');


VE = [(3^(1/2)*6^(1/2)*U^2)/12 + (3^(1/2)*6^(1/2)*(conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3))/2 - W/3 - (4*U)/3 + (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)/(2*conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3))) + (3^(1/2)*(conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3)) - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)/conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3)))*1i)/2 + 1)^2)/36 - (3^(1/2)*6^(1/2)*U)/4 - (3^(1/2)*6^(1/2)*W)/12 - (3^(1/2)*6^(1/2))/6 + (3^(1/2)*6^(1/2)*(4*U + W - 3)*(conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3))/2 - W/3 - (4*U)/3 + (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)/(2*conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3))) + (3^(1/2)*(conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3)) - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)/conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3)))*1i)/2 + 1))/36 + (3^(1/2)*6^(1/2)*U*W)/36, (6^(1/2)*U)/2 + (6^(1/2)*W)/6 + (6^(1/2)*(conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3))/2 - W/3 - (4*U)/3 + (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)/(2*conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3))) + (3^(1/2)*(conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3)) - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)/conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3)))*1i)/2 + 1))/6, 1]
alpha = VE(1);
beta = VE(2);
gamma = VE(3);
mag_alpha_sq = abs(alpha)^2;
mag_beta_sq = abs(beta)^2;
mag_gamma_sq = abs(gamma)^2;
sum = mag_alpha_sq + mag_beta_sq + mag_gamma_sq;
alpha1 = alpha/sqrt(sum);
beta1 = beta/sqrt(sum);
gamma1 = gamma/sqrt(sum);

E3 = (4*U)/3 + W/3 - (((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3)/2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)/(2*(((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3)) + (3^(1/2)*((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3) - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)/(((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3))*1i)/2 - 1
%ground state energy is third E value
%% Plot S1 vs U for dif U's
clc; clear all; close all;
J=1;
% Define the range of U values
U_values = 0:.1:100;

% Initialize arrays to store results
S1_values = zeros(size(U_values));
W_values = [0:1:10, 20:20:100, 1000];
%W_values = [1 2 4 8 16 32 64 128 256 512 1024];
% Initialize a figure for plotting
figure;
hold on;
Colorset = 0;
% Loop through each W value
for W = W_values
% Loop through each U value
for i = 1:length(U_values)
    U = U_values(i);
    
    % Recalculate alpha, beta, gamma for each U
    alpha = (2^(1/2)*U^2)/4 - (2^(1/2)*W)/4 - 2^(1/2)/2 - (3*2^(1/2)*U)/4 + (2^(1/2)*(conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3))/2 - W/3 - (4*U)/3 + (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)/(2*conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3))) + (3^(1/2)*(conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3)) - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)/conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3)))*1i)/2 + 1)^2)/12 + (2^(1/2)*U*W)/12 + (2^(1/2)*(4*U + W - 3)*(conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3))/2 - W/3 - (4*U)/3 + (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)/(2*conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3))) + (3^(1/2)*(conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3)) - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)/conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3)))*1i)/2 + 1))/12;
    beta = (6^(1/2)*U)/2 + (6^(1/2)*W)/6 + (6^(1/2)*(conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3))/2 - W/3 - (4*U)/3 + (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)/(2*conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3))) + (3^(1/2)*(conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3)) - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)/conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3)))*1i)/2 + 1))/6;
    gamma = 1;
 % Compute magnitudes squared
mag_alpha_sq = abs(alpha)^2;
mag_beta_sq = abs(beta)^2;
mag_gamma_sq = abs(gamma)^2;
sum = mag_alpha_sq + mag_beta_sq + mag_gamma_sq;
alpha1 = alpha/sqrt(sum);
beta1 = beta/sqrt(sum);
gamma1 = gamma/sqrt(sum);

lambda1 = (1/3)*abs(gamma1)^2;
lambda2 = (1/3)*abs(beta1)^2;
lambda3 = abs(alpha1)^2 + (1/3)*abs(beta1)^2;
lambda4 = (2/3)*abs(gamma1)^2 + (1/3)*abs(beta1)^2;
    % Calculate S1
    S1 = - lambda1*log(lambda1) - lambda2*log(lambda2)  - lambda3*log(lambda3) - lambda4*log(lambda4) ;
    % Store S1 value
    S1_values(i) = S1;

end

   % Calculate Colorset to transition from blue to red
    Colorset = Colorset + 1/(length(W_values)-1);
    
    if W == 0
        plot(log(U_values), S1_values, 'LineWidth', 4, 'Color', [0.1, 0.7, 1], 'DisplayName', 'U''/J = 0 (Unrestricted)');
    else
        % Calculate RGB components for the current color
        R = max(min(0.1 + 0.9*Colorset, 1), 0);  % Clamp R between 0 and 1
        G = max(min(0.7 - 0.5*Colorset, 1), 0);  % Clamp G between 0 and 1
        B = max(min(1 - Colorset, 1), 0);        % Clamp B between 0 and 1
        
        plot(log(U_values), S1_values, 'LineWidth', 1.5, 'Color', [R, G, B], ...
            'DisplayName', ['U''/J = ' num2str(W)],'LineStyle','-.');
    end
xlabel('log(U/J)','FontSize', 18);
ylabel('S1','FontSize', 18);
title('S1 Entanglement vs log(U/J)','FontSize', 18);
xlim([-2.3 4.57])
end

%% Graph for nmax = 2 case
%clc;clear all
syms J U W real
sqrt12 = sqrt(12);
sqrt6 = sqrt(6);
% 2x2 Hamiltonian Matrix
H = [0, -sqrt12;
    -sqrt12, -3 + U;];

% Compute eigenvalues and eigenvectors
[V, D] = eig(H);

% Extract eigenvalues
eigenvalues = diag(D);

% Display eigenvalues
disp('Eigenvalues:');
disp(eigenvalues);

% Display eigenvectors
disp('Eigenvectors:');
disp(V);

alphar = (3^(1/2)*((U^2 - 6*U + 57)^(1/2)/2 - U/2 + 3/2))/6 + (3^(1/2)*(U - 3))/6;
betar = 1;

%Graph
J=1; 
% Define the range of U values
U_values = 0:.1:100;
hold on;

% Loop through each U value
for i = 1:length(U_values)
    U = U_values(i);
    
    % Recalculate alpha, beta for each U
    alphar = (3^(1/2)*((U^2 - 6*U + 57)^(1/2)/2 - U/2 + 3/2))/6 + (3^(1/2)*(U - 3))/6;
    betar = 1;

 % Compute magnitudes squared
mag_alphar_sq = abs(alphar)^2;
mag_betar_sq = abs(betar)^2;
sum = mag_alphar_sq + mag_betar_sq;
alphar1 = alphar/sqrt(sum);
betar1 = betar/sqrt(sum);

lambda2 = (1/3)*abs(betar1)^2;
lambda3 = abs(alphar1)^2 + (1/3)*abs(betar1)^2;
lambda4 =  (1/3)*abs(betar1)^2;
    % Calculate S1
    S1R = - lambda2*log(lambda2)  - lambda3*log(lambda3) - lambda4*log(lambda4) ;
    % Store S1 value
    S1R_values(i) = S1R;


end

% Plot S1 vs U
plot(log(U_values), S1R_values, 'r--', 'LineWidth', 4, 'DisplayName', ['Restricted']);


% Particle Entanglement----------------------------------------------------------------------------------------
hold on;
J=1;
% Define the range of U values
U_values = 0:.1:100;

% Initialize arrays to store results
S1_values = zeros(size(U_values));
W_values = [0:1:10, 20:20:100, 1000];
%W_values = [1 2 4 8 16 32 64 128 256 512 1024];
% Initialize a figure for plotting
hold on;
Colorset = 0;
% Loop through each W value
for W = W_values
% Loop through each U value
for i = 1:length(U_values)
    U = U_values(i);
    
    % Recalculate alpha, beta, gamma for each U
    alpha = (2^(1/2)*U^2)/4 - (2^(1/2)*W)/4 - 2^(1/2)/2 - (3*2^(1/2)*U)/4 + (2^(1/2)*(conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3))/2 - W/3 - (4*U)/3 + (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)/(2*conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3))) + (3^(1/2)*(conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3)) - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)/conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3)))*1i)/2 + 1)^2)/12 + (2^(1/2)*U*W)/12 + (2^(1/2)*(4*U + W - 3)*(conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3))/2 - W/3 - (4*U)/3 + (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)/(2*conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3))) + (3^(1/2)*(conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3)) - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)/conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3)))*1i)/2 + 1))/12;
    beta = (6^(1/2)*U)/2 + (6^(1/2)*W)/6 + (6^(1/2)*(conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3))/2 - W/3 - (4*U)/3 + (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)/(2*conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3))) + (3^(1/2)*(conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3)) - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)/conj((((18*U + 6*W - ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 - (4*U + W - 3)^3/27)^2 - (3*U + W - (U*W)/3 - U^2 + (4*U + W - 3)^2/9 + 6)^3)^(1/2) - 6*W - 18*U + ((4*U + W - 3)*(9*U + 3*W - U*W - 3*U^2 + 18))/6 + (4*U + W - 3)^3/27)^(1/3)))*1i)/2 + 1))/6;
    gamma = 1;
 % Compute magnitudes squared
mag_alpha_sq = abs(alpha)^2;
mag_beta_sq = abs(beta)^2;
mag_gamma_sq = abs(gamma)^2;
sum = mag_alpha_sq + mag_beta_sq + mag_gamma_sq;
alpha1 = alpha/sqrt(sum);
beta1 = beta/sqrt(sum);
gamma1 = gamma/sqrt(sum);

X = ((2*gamma1*beta1)/(3*sqrt(6))) + (beta1^2)/6 + (2*alpha1*beta1)/(3*sqrt(3));
lambda1 = (1/3)+2*X;
lambda2 = (1/3)-X;
lambda3 = (1/3)-X;
    % Calculate S1
    S1 = - lambda1*log(lambda1) - lambda2*log(lambda2)  - lambda3*log(lambda3);
    % Store S1 value
    S1_values(i) = S1;

end

    % Calculate Colorset to transition from blue to red
    Colorset = Colorset + 1/(length(W_values)-1);
    
    if W == 0
        plot(log(U_values), S1_values, 'LineWidth', 4, 'Color', [0.1, 0.7, 1]);
    else
        % Calculate RGB components for the current color
        R = max(min(0.1 + 0.9*Colorset, 1), 0);  % Clamp R between 0 and 1
        G = max(min(0.7 - 0.5*Colorset, 1), 0);  % Clamp G between 0 and 1
        B = max(min(1 - Colorset, 1), 0);        % Clamp B between 0 and 1
        
        plot(log(U_values), S1_values, 'LineWidth', 1.5, 'Color', [R, G, B], 'LineStyle',':');
    end
end


%% Graph for nmax = 2 case
%clc;clear all
syms J U W real
alphar = (3^(1/2)*((U^2 - 6*U + 57)^(1/2)/2 - U/2 + 3/2))/6 + (3^(1/2)*(U - 3))/6;
betar = 1;

%Graph
J=1; 
% Define the range of U values
U_values = 0:.1:100;
hold on;

% Loop through each U value
for i = 1:length(U_values)
    U = U_values(i);
    
    % Recalculate alpha, beta for each U
    alphar = (3^(1/2)*((U^2 - 6*U + 57)^(1/2)/2 - U/2 + 3/2))/6 + (3^(1/2)*(U - 3))/6;
    betar = 1;

 % Compute magnitudes squared
mag_alphar_sq = abs(alphar)^2;
mag_betar_sq = abs(betar)^2;
sum = mag_alphar_sq + mag_betar_sq;
alphar1 = alphar/sqrt(sum);
betar1 = betar/sqrt(sum);
X = (betar1^2)/6 + (2*alphar1*betar1)/(3*sqrt(3));
lambda1 = (1/3)+2*X;
lambda2 = (1/3)-X;
lambda3 = (1/3)-X;
    % Calculate S1
    S1R = - lambda1*log(lambda1)  - lambda2*log(lambda2) - lambda3*log(lambda3) ;
    % Store S1 value
    S1R_values(i) = S1R;
end
% Plot S1 vs U
plot(log(U_values), S1R_values, 'r--', 'LineWidth', 4);

%removes particle from legend since same colors and stuff anyway
lines_with_name = findobj(gca, 'Type', 'line', '-not', 'DisplayName', '');
legend(lines_with_name,'Location','east','FontSize', 14);

% Add a text box for "Particle"
annotation('textbox', [0.15, 0.34, 0.05, 0.05], 'String', 'Particle Entanglement', 'FitBoxToText', 'on','FontSize', 20,'EdgeColor', 'none');

% Add a text box for "Spatial"
annotation('textbox', [0.15, 0.745, 0.05, 0.05], 'String', 'Spatial Entanglement', 'FitBoxToText', 'on','FontSize', 20,'EdgeColor', 'none');