%==========================================================================
% generate_divider_surface.m
%
% Function: Generate the reflector perturbation data file divider_surface.txt,
%           to be read by generate_reflector_stl.m and compute_H_from_surface.m.
%
% Description:
%   divider_surface.txt stores the perturbation values of 35 control nodes
%   (unit: m). These 35 nodes correspond to positions sampled every 3 points
%   from the 103-point full-resolution coordinates.
%   Parameter settings must be consistent with generate_reflector_stl.m.
%
% Perturbation modes (switched via the perturb_mode parameter):
%   'zero'  — no perturbation, all zeros (for verifying the standard paraboloid)
%   'sine'  — sinusoidal perturbation (for testing perturbation effect)
%
% Output file:
%   divider_surface.txt  — perturbation values of 35 control nodes
%                          (1 row x 35 columns, unit: m)
%
%==========================================================================

clc;
clear;

%% ========================================================================
%  Part 1: Parameter settings (must be consistent with generate_reflector_stl.m)
%% ========================================================================

% Reference scaling factor
scale = 3.8;

% Number of sample points per segment between control nodes
pts_per_segment = 3;

% Number of full-resolution sample points in x direction
n_x = 103;

% Total width in x direction (m)
x_span = (n_x - 1) * 3 / scale * 1e-3 / pts_per_segment;

% Number of control nodes (one every pts_per_segment points)
n_ctrl = (n_x - 1) / pts_per_segment + 1;   % = 35

% Full-resolution x coordinates (from positive to negative, consistent with generate_reflector_stl.m)
x_full = -linspace(-x_span/2, x_span/2, n_x);   % 1x103, unit: m

% x coordinates of control nodes (sampled every 3 points from full-resolution)
ctrl_node_idx = 1 : pts_per_segment : n_x;       % indices: 1,4,7,...,103
x_ctrl = x_full(ctrl_node_idx);                  % 1x35, unit: m

fprintf('Number of control nodes: %d\n', n_ctrl);
fprintf('x range: %.3f ~ %.3f mm\n', min(x_ctrl)*1e3, max(x_ctrl)*1e3);

%% ========================================================================
%  Part 2: Select perturbation mode
%% ========================================================================

% Toggle this parameter to select the perturbation type: 'zero' or 'sine'
perturb_mode = 'zero';

%% ========================================================================
%  Part 3: Generate perturbation data
%% ========================================================================

switch perturb_mode

    case 'zero'
        %------------------------------------------------------------------
        % No perturbation: all zeros
        % Purpose: verify standard paraboloid, used as the initial state for optimization
        %------------------------------------------------------------------
        surf_ctrl_pts = zeros(1, n_ctrl);
        fprintf('Perturbation mode: zero (no perturbation)\n');

    case 'sine'
        %------------------------------------------------------------------
        % Sinusoidal perturbation
        % Parameters:
        %   sine_amplitude  — perturbation amplitude (m), recommended not to exceed reflector depth limit
        %   sine_periods    — number of complete cycles within the reflector width
        %   sine_phase_deg  — initial phase (degrees)
        %------------------------------------------------------------------
        sine_amplitude = 0.1 / scale * 1e-3;   % perturbation amplitude (m), approx. 0.029 mm
        sine_periods   = 2;                      % number of sine periods
        sine_phase_deg = 0;                      % initial phase (degrees)

        surf_ctrl_pts = sine_amplitude ...
                        * sin(2 * pi * sine_periods * x_ctrl / x_span ...
                              + sine_phase_deg * pi / 180);

        fprintf('Perturbation mode: sine\n');
        fprintf('  Amplitude: %.4f um\n', sine_amplitude * 1e6);
        fprintf('  Periods: %d\n', sine_periods);

    otherwise
        error('Unknown perturbation mode: %s\nSupported modes: zero, sine', perturb_mode);
end

%% ========================================================================
%  Part 4: Write to file
%% ========================================================================

output_file = 'divider_surface.txt';
dlmwrite(output_file, surf_ctrl_pts, 'delimiter', ' ', 'precision', '%.10e');
fprintf('Perturbation file written to: %s\n', output_file);

%% ========================================================================
%  Part 5: Preview perturbation distribution
%% ========================================================================

figure(1);
stem(x_ctrl * 1e3, surf_ctrl_pts * 1e6, 'b', 'LineWidth', 1.5);
xlabel('x (mm)');
ylabel('Perturbation (um)');
title(sprintf('Reflector Control Node Perturbation (Mode: %s)', perturb_mode));
grid on;
