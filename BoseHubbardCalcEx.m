% Ryland Ross, June 2024
%------------------------------------- Spatial entanglement ---------------------------------------
clc; clear all; close all
syms U J

%% Finding e-values of B.H hamiltonian
A = [U -sqrt(2)*J 0; -sqrt(2)*J 0 -sqrt(2)*J; 0 -sqrt(2)*J U]; %Found A by hand, by writing out the matrix representation of the Hamiltonian (see pics in June 11th daily Log)

% Find the eigenvalues
eigenvalues = eig(A);
disp('The eigenvalues are:');
disp(eigenvalues);


%% Plotting the Spatial Entanglement
syms X
alpha = sqrt(4 / (X^2 + X * sqrt(X^2 + 16) + 16)); %found alpha/beta values by hand
gamma = alpha;
beta = sqrt(1 - (2 * alpha * alpha));
%this part is optional because they are all real values, but this is the more general case of using magnitude squared
alpha_conjugate = conj(alpha);
alpha_magnitude_squared = alpha * alpha_conjugate;
beta_conjugate = conj(beta);
beta_magnitude_squared = beta * beta_conjugate;
gamma_conjugate = conj(gamma);
gamma_magnitude_squared = gamma * gamma_conjugate;


% Define S1 using the magnitude squared of beta
S1_spacial = -((alpha_magnitude_squared * log(alpha_magnitude_squared)) + (beta_magnitude_squared * log(beta_magnitude_squared)) + (gamma_magnitude_squared * log(gamma_magnitude_squared)));
disp('The symbolic expression for S1 is:');
disp(S1_spacial);

% Define a range of U/J values
X_val = linspace(0, 10, 200); 

% Compute S1 for each U/J value
S1_spatial_values = zeros(size(X_val));
for i = 1:length(X_val)
    % Substitute the current U/J value into the expression for S1
    S1_spatial_values(i) = double(subs(S1_spacial, X, X_val(i))); 
end

% Plot S1 vs U/J
plot(X_val, S1_spatial_values,'LineWidth',2);
xlabel('U/J');
ylabel('S1 Entanglement');
title("Spatial and Particle Entanglement as a function of Interaction Strength");
grid on;
legend("Spatial");
hold on 

% ------------------------------------- Particle entanglement ---------------------------------------
lambda_plus = (1+2*sqrt(-4*alpha^4 + 2*alpha^2))/2;
lambda_minus = (1-2*sqrt(-4*alpha^4 + 2*alpha^2))/2;
S1_particle = -lambda_plus*log(lambda_plus)-lambda_minus*log(lambda_minus);

% Compute S1 for each U/J value
S1_particle_values = zeros(size(X_val));
for i = 1:length(X_val)
    % Substitute the current U/J value into the expression for S1
    alpha_val = double(subs(alpha, X, X_val(i)));
    lambda_plus_val = (1 + 2 * sqrt(-4 * alpha_val^4 + 2 * alpha_val^2)) / 2;
    lambda_minus_val = (1 - 2 * sqrt(-4 * alpha_val^4 + 2 * alpha_val^2)) / 2;
    
    % Ensure that we handle the log of zero properly
    if lambda_plus_val > 0
        term1 = -lambda_plus_val * log(lambda_plus_val);
    else
        term1 = 0;
    end
    
    if lambda_minus_val > 0
        term2 = -lambda_minus_val * log(lambda_minus_val);
    else
        term2 = 0;
    end
    S1_particle_values(i) = term1 + term2;
end

% Plot S1 vs U/J
plot(X_val, S1_particle_values,'LineWidth',2);
legend("Spatial", "Particle");