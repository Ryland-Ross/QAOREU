%Show the Relative error for L=N=3 for reduced ground state, vs non reduced
%as we increase U' and U
clc;
clear all;
close all

% Define constants
J = 1;
U_values = [0,3,5,10];
W = linspace(0, 1000, 1000); %W = U'

% Create a figure for relative error
figure;
hold on;

% Loop through each value of U
for i = 1:length(U_values)
    U = U_values(i);
    
    % Define the functions
    UnonZero = @(W) real((4*U)/3 - J + W/3 - (3^(1/2)*((J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)/(((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3) - (((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3))*1i)/2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)/(2*(((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3)) - (((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3)/2);
    E1 = @(W) real(U/2 - (3*J)/2 - (57*J^2 - 6*J*U + U^2)^(1/2)/2);
    
    % Compute corresponding values for UnonZero and E1
    real_UnonZero = arrayfun(@(w) UnonZero(w), W);
    real_E1 = arrayfun(@(w) E1(w), W);
    
    % Compute the absolute relative difference
    relative_diff = abs((real_UnonZero - real_E1) ./ real_E1);
    
    % Plot the relative difference
    plot(W, relative_diff, 'DisplayName', ['U = ' num2str(U)]);
end

xlabel('U''');
ylabel('Relative Error');
legend('show', 'Location', 'best');
title('Relative Difference between E_{gr} and E_g for Different U and U'' Values');
hold off;

%% Create a second figure for negative log of relative error
figure;
hold on;
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
end

xlabel('log(U''/J)');
ylabel('-log(absolute error)');
legend('show', 'Location', 'best');
title('Error between E_{restricted ground} and E_{ground}');
hold off;

%% plotting rel error vs 1/U'
figure;
hold on;
% Loop through each value of U
for i = 1:length(U_values)
    U = U_values(i);
    
    % Define the functions
    UnonZero = @(W) real((4*U)/3 - J + W/3 - (3^(1/2)*((J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)/(((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3) - (((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3))*1i)/2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)/(2*(((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3)) - (((18*J^2*U + 6*J^2*W - (4*U - 3*J + W)^3/27 - ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^2 - (J*W - (U*W)/3 + (4*U - 3*J + W)^2/9 + 6*J^2 - U^2 + 3*J*U)^3)^(1/2) - 18*J^2*U - 6*J^2*W + (4*U - 3*J + W)^3/27 + ((4*U - 3*J + W)*(18*J^2 + 9*J*U + 3*W*J - 3*U^2 - W*U))/6)^(1/3)/2);
    E1 = @(W) real(U/2 - (3*J)/2 - (57*J^2 - 6*J*U + U^2)^(1/2)/2);
    
    % Compute corresponding values for UnonZero and E1
    real_UnonZero = arrayfun(@(w) UnonZero(w), W);
    real_E1 = arrayfun(@(w) E1(w), W);
    
    % Compute the absolute relative difference
    relative_diff = abs((real_UnonZero - real_E1) ./ real_E1);
    
    % Plot the relative difference
    plot(1./W, relative_diff, 'DisplayName', ['U = ' num2str(U)]);
end

xlabel('1/U''');
xlim([0,0.04]);
ylabel('Relative Error');
legend('show', 'Location', 'best');
title('Relative Difference between E_{gr} and E_g for Different U and U'' Values');
hold off;
