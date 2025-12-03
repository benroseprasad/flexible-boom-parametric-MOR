clc; clear; 

% Import Matrices

fprintf('Loading matrices...\n');

K = spconvert(load('MatrixExport_STIF2.mtx'));
M = spconvert(load('MatrixExport_MASS2.mtx'));

N = size(K,1);



% Identify Base DOFs from abauqs penalty values

diagK = diag(K);
penaltyThreshold = 1e30;

fixedDOFs = find(diagK > penaltyThreshold);


% Interior DOFs
interiorDOFs = setdiff(1:N, fixedDOFs);

% Numerical Stabilization for eig

K(fixedDOFs,fixedDOFs) = 1e20;
M(fixedDOFs,fixedDOFs) = 1e-6;


% Constraint Modes


Kii = K(interiorDOFs,interiorDOFs);
Kib = K(interiorDOFs,fixedDOFs);

Psi = -Kii \ Kib;


% Compute Fixed-Interface Interior Modes

nInteriorModes = 100;


Mii = M(interiorDOFs,interiorDOFs);

opts.tol   = 1e-8;
opts.maxit = 500;

[V,~] = eigs(Kii, Mii, nInteriorModes, 'sm', opts);


% Craig–Bampton Transformation Matrix

nb = length(fixedDOFs);

T = [eye(nb),                  zeros(nb, nInteriorModes);
     Psi,                      V];


% Reduced Matrices

Kr = T' * K * T;
Mr = T' * M * T;




% Solve the Reduced Eigenproblem

[Vr_small, Dr_small] = eig(full(Kr), full(Mr));
freq_CB = sort( sqrt(diag(Dr_small)) / (2*pi) );


% Frequencies from Abaqus (Hz)

freq_abaqus = [ ...
    1.4139  1.8725  28.630  37.645  71.837 ...
    90.417 117.43 134.13 183.64 234.73 ...
    304.31 382.13 402.36 448.11 527.63 ...
    552.57 610.86 670.49 740.09 788.80 ];


% Abaqus vs Craig–Bampton Frequencies comparison



errors = zeros(size(freq_abaqus));
cb_match = zeros(size(freq_abaqus));

for i = 1:length(freq_abaqus)
    [~, idx] = min(abs(freq_CB - freq_abaqus(i)));
    cb_match(i) = freq_CB(idx);
    errors(i) = 100 * abs(freq_CB(idx) - freq_abaqus(i)) / freq_abaqus(i);

    fprintf('%2d\t%8.4f\t%8.4f\t%7.2f%%\n', ...
        i, freq_abaqus(i), freq_CB(idx), errors(i));
end




% Plot Modal Accuracy


% Frequency Error Plot
figure;
plot(1:length(errors), errors, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Mode Number', 'FontSize', 12);
ylabel('Error (%)', 'FontSize', 12);
title('CB vs Abaqus Frequency Error by Mode', 'FontSize', 14);
grid on; xlim([1 length(errors)]);

% Abaqus vs CB frequency
figure; hold on;
plot(1:length(freq_abaqus), freq_abaqus, 'o-', 'LineWidth', 2, ...
     'MarkerSize', 8, 'DisplayName', 'Abaqus');
plot(1:length(cb_match), cb_match, 's--', 'LineWidth', 2, ...
     'MarkerSize', 8, 'DisplayName', 'CB Reduced');

xlabel('Mode Number', 'FontSize', 12);
ylabel('Frequency (Hz)', 'FontSize', 12);
title('Abaqus vs Craig–Bampton Reduced Frequencies', 'FontSize', 14);
grid on; legend('Location','northwest');
xlim([1 length(freq_abaqus)]);

% Error percentage labels
for i = 1:length(errors)
    text(i, cb_match(i) + 0.05*max(cb_match), ...
         sprintf('%.1f%%', errors(i)), ...
         'HorizontalAlignment','center','FontSize',9);
end
