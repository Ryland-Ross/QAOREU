%% Find eigen values of the 3x3 matrix, found by simplifying 10x10 Hamiltonan for L=3=N
clc;
clear all;
syms U J W

H = [0, -sqrt(12)*J, 0;
    -sqrt(12)*J, -3*J + U, -sqrt(6)*J;
    0, -sqrt(6)*J, 3*U + W];

eigenvalues = eig(H)
%% Compare the above E-values to find lowest (ground state)
J = 100; U =400000; W =4000;
%E2 is the smallest

E1= (4*U)/3 - J + W/3 + (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)/(((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3) + (((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3)
E2= (4*U)/3 - J + W/3 - (3^(1/2)*((J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)/(((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3) - (((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3))*1i)/2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)/(2*(((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3)) - (((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3)/2
E3 = (4*U)/3 - J + W/3 + (3^(1/2)*((J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)/(((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3) - (((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3))*1i)/2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)/(2*(((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3)) - (((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3)/2


%% Same thing, but now make U'=0
clear all;clc;
syms U J
H = [0, -sqrt(12)*J, 0;
    -sqrt(12)*J, -3*J + U, -sqrt(6)*J;
    0, -sqrt(6)*J, 3*U ];

eigenvalues = eig(H)
%% Compare the above E-values to find lowest (ground state)
% e2 is smallest
J = 100; U =40000;
E1 =                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   (4*U)/3 - J + ((3*J - 4*U)^2/9 + 6*J^2 - U^2 + 3*J*U)/((((3*J - 4*U)^3/27 + 18*J^2*U + ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6)^2 - ((3*J - 4*U)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6 - (3*J - 4*U)^3/27)^(1/3) + ((((3*J - 4*U)^3/27 + 18*J^2*U + ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6)^2 - ((3*J - 4*U)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6 - (3*J - 4*U)^3/27)^(1/3)
E2 =(4*U)/3 - J - ((3*J - 4*U)^2/9 + 6*J^2 - U^2 + 3*J*U)/(2*((((3*J - 4*U)^3/27 + 18*J^2*U + ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6)^2 - ((3*J - 4*U)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6 - (3*J - 4*U)^3/27)^(1/3)) - (3^(1/2)*(((3*J - 4*U)^2/9 + 6*J^2 - U^2 + 3*J*U)/((((3*J - 4*U)^3/27 + 18*J^2*U + ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6)^2 - ((3*J - 4*U)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6 - (3*J - 4*U)^3/27)^(1/3) - ((((3*J - 4*U)^3/27 + 18*J^2*U + ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6)^2 - ((3*J - 4*U)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6 - (3*J - 4*U)^3/27)^(1/3))*1i)/2 - ((((3*J - 4*U)^3/27 + 18*J^2*U + ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6)^2 - ((3*J - 4*U)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6 - (3*J - 4*U)^3/27)^(1/3)/2
E3 = (4*U)/3 - J - ((3*J - 4*U)^2/9 + 6*J^2 - U^2 + 3*J*U)/(2*((((3*J - 4*U)^3/27 + 18*J^2*U + ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6)^2 - ((3*J - 4*U)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6 - (3*J - 4*U)^3/27)^(1/3)) + (3^(1/2)*(((3*J - 4*U)^2/9 + 6*J^2 - U^2 + 3*J*U)/((((3*J - 4*U)^3/27 + 18*J^2*U + ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6)^2 - ((3*J - 4*U)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6 - (3*J - 4*U)^3/27)^(1/3) - ((((3*J - 4*U)^3/27 + 18*J^2*U + ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6)^2 - ((3*J - 4*U)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6 - (3*J - 4*U)^3/27)^(1/3))*1i)/2 - ((((3*J - 4*U)^3/27 + 18*J^2*U + ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6)^2 - ((3*J - 4*U)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6 - (3*J - 4*U)^3/27)^(1/3)/2

%% Finding e-values for 2x2 matrix, corresponding to nmax = 2 for L=3=N
clc;
clear all;
syms U J

H = [0, -sqrt(12)*J;
    -sqrt(12)*J, -3*J+U];

eigenvalues = eig(H)
%E1 is the smallest
J = 1; U =4; 
E1 = U/2 - (3*J)/2 - (57*J^2 - 6*J*U + U^2)^(1/2)/2
E2 = U/2 - (3*J)/2 + (57*J^2 - 6*J*U + U^2)^(1/2)/2
%% Compare nmax = 2 evalues graphically
clc; clear all; close all;
% Define constants
U = linspace(0, 100, 100);  % Range of U values
J = linspace(0, 100, 100);  % Range of J values

% Calculate E1 and E2
E1 = @(U, J) real(U/2 - (3*J)/2 - sqrt(57*J^2 - 6*J*U + U^2)/2);
E2 = @(U, J) real(U/2 - (3*J)/2 + sqrt(57*J^2 - 6*J*U + U^2)/2);

% Create a grid of U and J values
[U_grid, J_grid] = meshgrid(U, J);

% Calculate E1 and E2 over the grid
E1_grid = E1(U_grid, J_grid);
E2_grid = E2(U_grid, J_grid);

% Plotting
figure;

% Contour plot for E1
subplot(1, 2, 1);
contourf(U_grid, J_grid, E1_grid, 20, 'LineWidth', 1.5);
colorbar;
xlabel('U');
ylabel('J');
title('Contour Plot of E1');
grid on;

% Contour plot for E2
subplot(1, 2, 2);
contourf(U_grid, J_grid, E2_grid, 20, 'LineWidth', 1.5);
colorbar;
xlabel('U');
ylabel('J');
title('Contour Plot of E2');
grid on;
sgtitle('L=3=N, nmax = 2');


%% Graphing the groud states, increasing U'
clc;
clear all;
close all

% Define constants
J = 1;

% Define the values of U
U_values = [0, 1, 3, 5, 7, 8, 9, 10];

% Generate values for W
W = linspace(0, 2500, 1000);

% Create a figure for subplots
figure;

% Loop through each value of U
for i = 1:length(U_values)
    U = U_values(i);
    
    % Define the functions
    UnonZero = @(W) real((4*U)/3 - J + W/3 - (3^(1/2)*((J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)/(((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3) - (((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3))*1i)/2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)/(2*(((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3)) - (((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3)/2);
    
    Uzero = @(W) real((4*U)/3 - J - ((3*J - 4*U)^2/9 + 6*J^2 - U^2 + 3*J*U)/(2*((((3*J - 4*U)^3/27 + 18*J^2*U + ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6)^2 - ((3*J - 4*U)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6 - (3*J - 4*U)^3/27)^(1/3)) - (3^(1/2)*(((3*J - 4*U)^2/9 + 6*J^2 - U^2 + 3*J*U)/((((3*J - 4*U)^3/27 + 18*J^2*U + ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6)^2 - ((3*J - 4*U)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6 - (3*J - 4*U)^3/27)^(1/3) - ((((3*J - 4*U)^3/27 + 18*J^2*U + ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6)^2 - ((3*J - 4*U)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6 - (3*J - 4*U)^3/27)^(1/3))*1i)/2 - ((((3*J - 4*U)^3/27 + 18*J^2*U + ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6)^2 - ((3*J - 4*U)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6 - (3*J - 4*U)^3/27)^(1/3)/2);
    
    E1 = @(W) real(U/2 - (3*J)/2 - (57*J^2 - 6*J*U + U^2)^(1/2)/2);
    
    % Compute corresponding values for UnonZero, Uzero, and E1
    real_UnonZero = arrayfun(@(w) UnonZero(w), W);
    real_Uzero = arrayfun(@(w) Uzero(w), W);
    real_E1 = arrayfun(@(w) E1(w), W);
    
    % Plot in a subplot
    subplot(4, 2, i);
    plot(W, real_UnonZero, 'b', W, real_Uzero, 'r', W, real_E1, 'm');
    xlabel('U''');
ylabel('Ground State Energy');
legend('3x3, increasing U'' (W)', '3x3, U'' = 0', '2x2', 'Location', 'best');
grid on;
    title(['U = ' num2str(U)]);
end

% Add a main title for the entire figure
sgtitle('Ground State Energy of 2x2 vs 3x3 Hamiltonian against U'', for Different U values. J=1');

%% Plot just for U = 10
clc;
clear all;
close all

% Define constants
J = 1;

% Define the values of U
U_values = [10];

% Generate values for W
W = linspace(0, 2500, 1000);


% Loop through each value of U
for i = 1:length(U_values)
    U = U_values(i);
    
    % Define the functions
    UnonZero = @(W) real((4*U)/3 - J + W/3 - (3^(1/2)*((J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)/(((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3) - (((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3))*1i)/2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)/(2*(((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3)) - (((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3)/2);
    
    Uzero = @(W) real((4*U)/3 - J - ((3*J - 4*U)^2/9 + 6*J^2 - U^2 + 3*J*U)/(2*((((3*J - 4*U)^3/27 + 18*J^2*U + ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6)^2 - ((3*J - 4*U)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6 - (3*J - 4*U)^3/27)^(1/3)) - (3^(1/2)*(((3*J - 4*U)^2/9 + 6*J^2 - U^2 + 3*J*U)/((((3*J - 4*U)^3/27 + 18*J^2*U + ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6)^2 - ((3*J - 4*U)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6 - (3*J - 4*U)^3/27)^(1/3) - ((((3*J - 4*U)^3/27 + 18*J^2*U + ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6)^2 - ((3*J - 4*U)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6 - (3*J - 4*U)^3/27)^(1/3))*1i)/2 - ((((3*J - 4*U)^3/27 + 18*J^2*U + ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6)^2 - ((3*J - 4*U)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6 - (3*J - 4*U)^3/27)^(1/3)/2);
    
    E1 = @(W) real(U/2 - (3*J)/2 - (57*J^2 - 6*J*U + U^2)^(1/2)/2);
    
    % Compute corresponding values for UnonZero, Uzero, and E1
    real_UnonZero = arrayfun(@(w) UnonZero(w), W);
    real_Uzero = arrayfun(@(w) Uzero(w), W);
    real_E1 = arrayfun(@(w) E1(w), W);
    
 
figure;
plot(W, real_UnonZero, 'Color', [0.6, 0.2, .8], 'LineWidth', 2);
hold on;
plot(W, real_Uzero, 'Color', [0.1, 0.7, 1], 'LineWidth', 2);
plot(W, real_E1, 'r', 'LineWidth', 2);
hold off;
    xlabel('U''/J');
ylabel('Ground State Energy');
legend('Unrestricted, increasing U''', 'Unrestricted, U'' = 0', 'Restricted', 'Position', [0.68, 0.688, 0.1, 0.1]);
end

% Add a main title for the entire figure
title('Ground State Energy, U=10, J=1');
