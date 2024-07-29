% Define constants
clc; close all;
J = 1;
U_values = [10];

% Generate values for W
W = linspace(0, 2500, 1000);
    
    % Define the functions for the first plot
    UnonZero = @(W) real((4*U)/3 - J + W/3 - (3^(1/2)*((J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)/(((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3) - (((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3))*1i)/2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)/(2*(((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3)) - (((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3)/2);
    Uzero = @(W) real((4*U)/3 - J - ((3*J - 4*U)^2/9 + 6*J^2 - U^2 + 3*J*U)/(2*((((3*J - 4*U)^3/27 + 18*J^2*U + ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6)^2 - ((3*J - 4*U)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6 - (3*J - 4*U)^3/27)^(1/3)) - (3^(1/2)*(((3*J - 4*U)^2/9 + 6*J^2 - U^2 + 3*J*U)/((((3*J - 4*U)^3/27 + 18*J^2*U + ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6)^2 - ((3*J - 4*U)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6 - (3*J - 4*U)^3/27)^(1/3) - ((((3*J - 4*U)^3/27 + 18*J^2*U + ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6)^2 - ((3*J - 4*U)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6 - (3*J - 4*U)^3/27)^(1/3))*1i)/2 - ((((3*J - 4*U)^3/27 + 18*J^2*U + ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6)^2 - ((3*J - 4*U)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - ((3*J - 4*U)*(18*J^2 + 9*J*U - 3*U^2))/6 - (3*J - 4*U)^3/27)^(1/3)/2);
    E1 = @(W) real(U/2 - (3*J)/2 - (57*J^2 - 6*J*U + U^2)^(1/2)/2);
    
    % Compute corresponding values for UnonZero, Uzero, and E1
    real_UnonZero = arrayfun(@(w) UnonZero(w), W);
    real_Uzero = arrayfun(@(w) Uzero(w), W);
    real_E1 = arrayfun(@(w) E1(w), W);
    
    % Plot the first graph
    figure(1);
    plot(W, real_UnonZero, 'Color', [0.6, 0.2, .8], 'LineWidth', 2);
    hold on;
    plot(W, real_Uzero, 'Color', [0.1, 0.7, 1], 'LineWidth', 2);
    plot(W, real_E1, 'r', 'LineWidth', 2);
    hold off;
    xlabel('U''/J');
    ylabel('Ground State Energy');ylim([-1.455 -1.4235])
    legend('Unrestricted, increasing U''', 'Unrestricted, U'' = 0', 'Restricted', 'Position', [0.68, 0.76, 0.1, 0.1]);
    title(['Ground State Energy, U=' num2str(U) ', J=' num2str(J)]);
        inset_position = [0.22 0.27 0.67 0.42]; % [left bottom width height]
    axes('Position', inset_position);
   U_values = [0,3,5,10];
W = linspace(0, 1000, 1000); %W = U'
% Loop through each value of U
for i = 1:length(U_values)
    U = U_values(i);
    
    % Recalculate relative difference for each U
    UnonZero = @(W) real((4*U)/3 - J + W/3 - (3^(1/2)*((J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)/(((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3) - (((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3))*1i)/2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)/(2*(((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3)) - (((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3)/2);
    E1 = @(W) real(U/2 - (3*J)/2 - (57*J^2 - 6*J*U + U^2)^(1/2)/2);
    
    % Compute corresponding values for UnonZero and E1
    real_UnonZero = arrayfun(@(w) UnonZero(w), W);
    real_E1 = arrayfun(@(w) E1(w), W);
    
    % Compute the absolute relative difference
    relative_diff = abs((real_UnonZero - real_E1) ./ real_E1);
    
    % Plot the negative log of relative difference
    plot(log(W), -log(relative_diff), 'DisplayName', ['U = ' num2str(U)],'LineWidth',2);
    hold on;
end
legend('show', 'Location', 'best');
    xlabel('log(U''/J)');
    ylabel('-log(absolute relative error)');
    title('Error between E_{restricted ground} and E_{ground}');
