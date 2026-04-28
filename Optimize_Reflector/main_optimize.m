%==========================================================================
% main_optimize.m  —— Top-level reflector surface optimization script
%
% Function: Read the initial reflector surface, call the evolutionary
%           algorithm to optimize the surface perturbation values, display
%           the convergence curve in real time, and output the final
%           magnetic field distribution after optimization.
%
% Call hierarchy:
%   main_optimize.m
%     └── optimize_reflector_surface.m  (evolutionary algorithm)
%           └── calc_H_field.m          (electromagnetic field computation)
%           └── eval_field_correlation.m        (fitness evaluation)
%
% Dependency files (must be in the same directory as this script):
%   divider_surface.txt  — initial reflector surface (35 control nodes; all zeros = standard paraboloid)
%   Divider_C.txt        — target field pattern (1x153 complex vector)
%   Weight.txt           — weighting window function (1x153, 1 in center region)
%   Divider_mag.txt      — reference magnitude (optional, for comparison plot)
%   Divider_phase.txt    — reference phase (optional, for comparison plot)
%
% Output:
%   divider_surface.txt  — optimized best reflector surface (overwritten)
%   Figure 1: real-time convergence curve (error vs. generation)
%   Figure 2: final H3-plane H-field magnitude (dB) and phase distribution
%==========================================================================

clc;
clear;
tic;

%% ========================================================================
%  Part 1: Read the initial reflector surface
%% ========================================================================

init_file = 'divider_surface.txt';
if ~isfile(init_file)
    error('Initial reflector surface file not found: %s\nPlease run generate_divider_surface.m first to generate it.', init_file);
end

init_surf = dlmread(init_file);
fprintf('Initial reflector surface loaded: %s (nodes: %d)\n', init_file, numel(init_surf));

%% ========================================================================
%  Part 2: Set optimization parameters
%% ========================================================================

params.n_generations = 300;              % number of generations
params.pop_size      = 10;              % population size
params.surf_depth    = 1 / 3.8 * 1e-3; % perturbation depth limit (m)
params.Fe            = 0.1;             % DE crossover weight
params.Fm            = 0.1;             % mutation strength
params.mut_num       = 2;               % number of nodes mutated per step

%% ========================================================================
%  Part 3: Create the convergence curve figure window and pass it to the
%          optimization function
%% ========================================================================

fig_conv = figure(1);
clf(fig_conv);
xlabel('Generation');
ylabel('Error (100 - Fitness)');
title('Reflector Surface Optimization Convergence (running...)');
grid on;
drawnow;

params.fig_handle = fig_conv;   % pass figure handle to the optimization function for real-time updates

%% ========================================================================
%  Part 4: Run optimization
%% ========================================================================

fprintf('\nStarting optimization: %d generations, population size %d...\n\n', params.n_generations, params.pop_size);

[best_surf, best_fitness, fitness_history] = ...
    optimize_reflector_surface(init_surf, params);

fprintf('\nOptimization complete. Best fitness: %.4f\n', best_fitness);

%% ========================================================================
%  Part 5: Save the best reflector surface
%% ========================================================================

dlmwrite('divider_surface.txt', best_surf, 'delimiter', ' ', 'precision', '%.10e');
fprintf('Best reflector surface saved to: divider_surface.txt\n');

%% ========================================================================
%  Part 6: Update the convergence curve title
%% ========================================================================

figure(fig_conv);
title(sprintf('Reflector Surface Optimization Convergence (done, best fitness: %.4f)', best_fitness));

%% ========================================================================
%  Part 7: Compute and plot the final H3-plane magnetic field distribution
%% ========================================================================

H3_final = calc_H_field(best_surf);

% Magnitude (dB, clipped to -30 dB)
H3_abs  = abs(H3_final);
H3_max  = max(H3_abs);
H3_dB   = 20 * log10(H3_abs / H3_max);
H3_dB   = max(H3_dB, -30);

% Phase (degrees, with 100 deg reference phase offset, consistent with original program)
Phase_offset = 100;
H3_phase = angle(H3_final * exp(sqrt(-1) * Phase_offset * pi / 180)) * 180 / pi;

% Plot
figure(2);
clf;

subplot(2, 1, 1);
plot(H3_dB, 'b', 'LineWidth', 2);
hold on;
if isfile('Divider_mag.txt')
    ref_mag = dlmread('Divider_mag.txt');
    plot(ref_mag, 'r--', 'LineWidth', 1.5);
    legend('Optimized', 'Reference');
end
xlabel('H3 Plane Sample Index');
ylabel('Magnitude (dB)');
title('H3 Plane Magnetic Field Magnitude Distribution');
ylim([-30, 0]);
grid on;

subplot(2, 1, 2);
plot(H3_phase, 'b', 'LineWidth', 2);
hold on;
if isfile('Divider_phase.txt')
    ref_phase = dlmread('Divider_phase.txt');
    plot(ref_phase, 'r--', 'LineWidth', 1.5);
    legend('Optimized', 'Reference');
end
xlabel('H3 Plane Sample Index');
ylabel('Phase (deg)');
title('H3 Plane Magnetic Field Phase Distribution');
grid on;

%% ========================================================================
%  Part 8: Print final statistics
%% ========================================================================

n3x = 153;
scale = 3.8;
pts_per = 3;
x3max = (n3x - 1) * 3 / scale * 1e-3 / pts_per;
d3x   = x3max / (n3x - 1);

P_output = 0.5 * 377 * sum(abs(H3_final).^2) * d3x;
fprintf('Output power (normalized): %.6f\n', P_output);

elapsed = toc;
fprintf('Total elapsed time: %.1f s\n', elapsed);
