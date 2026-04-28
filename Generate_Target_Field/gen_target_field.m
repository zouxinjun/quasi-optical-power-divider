%==========================================================================
% gen_target_field.m
%
% Function: Generate the target magnetic field distribution of a power
%           divider based on 8 beam amplitude weights, and save the results
%           as data files for use by the optimization program.
%
% Principle:
%   1. Compute the magnetic field distribution at a single horn aperture (TE10 mode)
%   2. Concatenate the fields of 8 horns with weights to form the 8-way
%      power divider target pattern
%   3. Embed the target pattern into the center region of the H3 observation
%      plane (153 points)
%   4. Save the complex field, magnitude (dB), and phase to files
%
% Output files:
%   Divider_C.txt      — target field pattern (1x153 complex vector)
%   Divider_mag.txt    — target field magnitude (1x153, dB, clipped to -30 dB)
%   Divider_phase.txt  — target field phase (1x153, degrees)
%
%==========================================================================

clc;
clear;
tic;

%% ========================================================================
%  Part 1: 8-beam amplitude weight settings
%  Each weight defines the output amplitude ratio of one horn antenna.
%  Adjust manually as needed.
%% ========================================================================

% 8-beam amplitude weights (dimensionless, relative ratios)
% Left to right correspond to output ports 1 through 8
beam_weights = [1.25, 1.25, 0.7, 0.8, 0.5, 1.1, 1.55, 1.7];
%beam_weights = [1, 1, 1, 1, 1, 1, 1, 1];

n_beams = numel(beam_weights);   % number of beams = 8

%% ========================================================================
%  Part 2: Basic electromagnetic parameters
%% ========================================================================

freq_GHz = 380;
lambda   = 300 / freq_GHz * 1e-3;   % free-space wavelength (m)
beta     = 2 * pi / lambda;          % phase constant (rad/m)
j_unit   = sqrt(-1);

%% ========================================================================
%  Part 3: Single horn aperture parameters
%% ========================================================================

pts_per_lambda = 3;   % number of sample points per wavelength (sampling density)

% Number of sample points at the horn aperture (x: E-plane wide wall, y: H-plane narrow wall)
n_horn_x = 7;    % wide-wall sample points
n_horn_y = 12;   % narrow-wall sample points (only x direction is used; noted here for reference)

% Physical dimensions of the horn aperture (derived from sample count and density)
horn_width  = (n_horn_x - 1) * lambda / pts_per_lambda;   % wide-wall aperture width (m)
horn_height = (n_horn_y - 1) * lambda / pts_per_lambda;   % narrow-wall aperture height (m)

% Horn aperture coordinates (uniformly spaced, centered at origin)
x_horn = linspace(-horn_width/2, horn_width/2, n_horn_x);   % wide-wall coordinates (m)

% Horn geometry parameters
horn_length = 9 / 3.8 * 1e-3;   % horn aperture length (m)
wg_width    = 0.508e-3;          % waveguide wide-wall dimension (m)
wg_height   = 0.254e-3;          % waveguide narrow-wall dimension (m)

% Equivalent phase-center radii of curvature (for computing aperture phase distribution)
R_horn_x = horn_length / (1 - (wg_width  / horn_width));    % wide-wall equivalent radius (m)
R_horn_y = horn_length / (1 - (wg_height / horn_height));   % narrow-wall equivalent radius (m, reserved)

E0 = 377;   % initial electric field amplitude (V/m, normalized to free-space wave impedance)

%% ========================================================================
%  Part 4: Compute the magnetic field distribution at a single horn
%          aperture (TE10 mode)
%
%  Model: Hx = -(E0/377) * cos(pi*x/DH) * exp(-j*beta*x^2/(2*R))
%  Only the x direction (wide wall) is considered.
%% ========================================================================

H_horn = zeros(1, n_horn_x);
for j = 1 : n_horn_x
    x = x_horn(j);
    % Cosine amplitude distribution (TE10 mode) x spherical-wave phase (equivalent phase center)
    H_horn(j) = -E0 / 377 * cos(pi * x / horn_width) ...
                * exp(-j_unit * beta * x^2 / (2 * R_horn_x));
end

% Apply reference phase offset (-180 deg) to match the original program
H_horn = H_horn * exp(j_unit * (-180) * pi / 180);

% Magnitude (dB) and phase (degrees) of the single horn aperture field, for plotting reference
H_horn_max   = max(abs(H_horn));
H_horn_dB    = 20 * log10(abs(H_horn) / H_horn_max);
H_horn_phase = angle(H_horn) * 180 / pi;

% Clip edge points (edge values are negligible; avoids log computation anomalies)
H_horn_dB(1)         = -30;
H_horn_dB(n_horn_x)  = -30;

%% ========================================================================
%  Part 5: Concatenate 8 beams to build the target field pattern
%
%  Concatenation rule: adjacent horns share their boundary sample point.
%  Total points = n_horn_x * n_beams - (n_beams - 1)
%  Each beam field is multiplied by the corresponding amplitude weight.
%% ========================================================================

n_total = n_horn_x * n_beams - (n_beams - 1);   % total sample points after concatenation = 49

H_target_raw = zeros(1, n_total);   % concatenated target field pattern (complex)

for b = 1 : n_beams
    % Starting index of beam b in the concatenated array
    start_idx = (b - 1) * (n_horn_x - 1) + 1;
    for j = 1 : n_horn_x
        H_target_raw(start_idx + j - 1) = H_horn(j) * beam_weights(b);
    end
end

% Magnitude (dB) and phase (degrees) of the concatenated field
H_raw_max   = max(abs(H_target_raw));
H_raw_dB    = 20 * log10(abs(H_target_raw) / H_raw_max);
H_raw_dB    = max(H_raw_dB, -30);
H_raw_phase = angle(H_target_raw) * 180 / pi;

%% ========================================================================
%  Part 6: Embed the target field pattern into the center of the H3
%          observation plane (153 points)
%
%  The H3 observation plane has 153 sample points. The target pattern
%  (49 points) is placed at the center; both ends are zero-padded,
%  indicating no target field requirement at the observation plane edges.
%% ========================================================================

n3x = 153;   % total sample points on H3 observation plane (consistent with calc_H_field.m)

H_target = zeros(1, n3x);   % initialized to all zeros (no target field at edges)

% Compute the starting index for center-aligned embedding
offset = (n3x - n_total) / 2;   % number of zero-padded points on the left
for j = 1 : n_total
    H_target(offset + j) = H_target_raw(j);
end

% Magnitude (dB) and phase (degrees) of the embedded target field
H_target_max   = max(abs(H_target));
H_target_dB    = 20 * log10(abs(H_target) / H_target_max);
H_target_dB    = max(H_target_dB, -30);
H_target_phase = angle(H_target) * 180 / pi;

%% ========================================================================
%  Part 7: Preview plot
%% ========================================================================

figure(1);

subplot(2, 3, 1);
plot(x_horn * 1e3, H_horn_dB, 'b', 'LineWidth', 1.5);
xlabel('x (mm)');
ylabel('Magnitude (dB)');
title('Single Horn Aperture Field Magnitude');
grid on;

subplot(2, 3, 4);
plot(x_horn * 1e3, H_horn_phase, 'r', 'LineWidth', 1.5);
xlabel('x (mm)');
ylabel('Phase (deg)');
title('Single Horn Aperture Field Phase');
grid on;

subplot(2, 3, 2);
plot(H_raw_dB, 'b', 'LineWidth', 1.5);
xlabel('Sample Index');
ylabel('Magnitude (dB)');
title('8-Beam Concatenated Target Field Magnitude');
grid on;

subplot(2, 3, 5);
plot(H_raw_phase, 'r', 'LineWidth', 1.5);
xlabel('Sample Index');
ylabel('Phase (deg)');
title('8-Beam Concatenated Target Field Phase');
grid on;

subplot(2, 3, 3);
plot(H_target_dB, 'b', 'LineWidth', 1.5);
xlabel('H3 Plane Sample Index');
ylabel('Magnitude (dB)');
title('H3 Plane Target Field Magnitude (153 points)');
grid on;

subplot(2, 3, 6);
plot(H_target_phase, 'r', 'LineWidth', 1.5);
xlabel('H3 Plane Sample Index');
ylabel('Phase (deg)');
title('H3 Plane Target Field Phase (153 points)');
grid on;

%% ========================================================================
%  Part 8: Save data files
%% ========================================================================

dlmwrite('Divider_C.txt',     H_target,       'delimiter', ' ', 'precision', 4);
dlmwrite('Divider_mag.txt',   H_target_dB,    'delimiter', ' ', 'precision', 4);
dlmwrite('Divider_phase.txt', H_target_phase, 'delimiter', ' ', 'precision', 4);

fprintf('Target field pattern saved:\n');
fprintf('  Divider_C.txt      — complex field (1x%d)\n', n3x);
fprintf('  Divider_mag.txt    — magnitude (dB)\n');
fprintf('  Divider_phase.txt  — phase (degrees)\n');

toc;
