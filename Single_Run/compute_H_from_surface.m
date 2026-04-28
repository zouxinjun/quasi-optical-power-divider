%==========================================================================
% compute_H_from_surface.m
%
% Function: Read the reflector surface file divider_surface.txt, compute
%           the magnetic field (H-field) distribution of a terahertz-band
%           power divider using the Physical Optics (PO) method, and output
%           the results.
%
% Principle:
%   1. The horn antenna aperture (H1 plane) generates the initial H-field
%   2. The field at H1 propagates through free space and illuminates the
%      parabolic + perturbed reflector surface (H2 plane)
%   3. H2 acts as a secondary radiation source and propagates again to the
%      observation plane (H3 plane)
%   4. Output the H-field magnitude (dB) and phase distribution at H3
%
% Input file:
%   divider_surface.txt  — reflector surface perturbation values (35 nodes, 1 row)
%
% Output:
%   Console printout of power, fitness, and other metrics
%   Plot of H3-plane magnitude (dB) and phase
%   Divider_surface_backup.txt — backup of the loaded surface data
%
%==========================================================================

clc;
clear;
tic;

%% ========================================================================
%  Part 1: Basic physical parameters
%% ========================================================================

freq_GHz   = 380;                        % operating frequency (GHz)
lambda     = 300 / freq_GHz * 1e-3;     % free-space wavelength (m)
beta       = 2 * pi / lambda;            % phase constant (rad/m)
j_unit     = sqrt(-1);                   % imaginary unit

% Reference scaling factor (original design based on 3.8 GHz, scaled proportionally)
scale = 3.8;

%% ========================================================================
%  Part 2: Geometric parameters
%% ========================================================================

% --- Maximum reflector perturbation depth ---
surf_depth = 1 / scale * 1e-3;          % maximum perturbation depth (m)

% --- Incidence angle (between horn axis and reflector normal) ---
theta      = 45 * pi / 180;             % 45 deg

% Precomputed trigonometric values
cos_t  = cos(theta);
sin_t  = sin(theta);
cos_2t = cos(2 * theta);
sin_2t = sin(2 * theta);

% --- Paraboloid focal length ---
focus = 0.11 / scale;                   % focal length (m)

% --- Coordinate translations ---
z2_shift = 66 / 3.6 * 1e-3;            % H2-plane z-axis translation (m)
z_HFSS   = 10 / scale * 1e-3;          % HFSS simulation reference plane offset (m)
z3_shift = (z2_shift * sin_2t * tan(2*theta) - z_HFSS) * cos_2t;
x3_shift = 50 / scale * 1e-3;          % H3-plane x-axis translation (m)

%% ========================================================================
%  Part 3: Horn antenna aperture (H1 plane) parameters
%% ========================================================================

horn_length = 22 / scale * 1e-3;       % horn aperture length (m)
wg_width    = 0.508e-3;                 % waveguide wide-wall dimension (m)
wg_height   = 0.254e-3;                 % waveguide narrow-wall dimension (m)
horn_aperture = 20 / scale * 1e-3;     % horn aperture width (m)

% Equivalent phase-center radius of the horn (for aperture phase distribution)
R_horn = horn_length / (1 - (wg_width / horn_aperture));

% H1-plane sampling
n1x   = 15;                             % number of sample points on H1 plane
x1max = horn_aperture / 2;             % half-width of H1 plane
dx1   = horn_aperture / (n1x - 1);     % sample spacing on H1 plane
ds1   = dx1;                            % integration element length on H1 plane
x1_coords = linspace(-x1max, x1max, n1x);  % x coordinates on H1 plane

E0 = 377;  % initial electric field amplitude (normalized to free-space wave impedance 377 ohm)

%% ========================================================================
%  Part 4: Compute the H-field distribution at H1 plane
%  Model: TE10-mode horn aperture field with spherical-wave phase factor
%% ========================================================================

H1 = zeros(1, n1x);
for idx = 1 : n1x
    x = x1_coords(idx);
    % Amplitude: cosine distribution (TE10 mode); phase: spherical-wave phase
    H1(idx) = -E0 / 377 * cos(pi * x / horn_aperture) ...
              * exp(-j_unit * beta * x^2 / (2 * R_horn));
end

% Print normalized input power at H1 plane
P_input = 0.5 * sum(abs(H1).^2) * ds1;
fprintf('H1 plane normalized input power P_input = %.6f\n', P_input);

%% ========================================================================
%  Part 5: H2 plane (reflector surface) parameters
%% ========================================================================

n2x   = 103;                                        % number of sample points on H2 plane
x2max = (n2x - 1) * 3 / scale * 1e-3 / 3;         % half-width of H2 plane (m)
d2x   = x2max / (n2x - 1);                         % sample spacing on H2 plane
ds2   = d2x;                                        % integration element length on H2 plane
n2x_half = (n2x - 1) / 2;                          % half number of points on H2 plane

% H2-plane x coordinates (sign-flipped to match coordinate system)
x2_coords = -linspace(-x2max/2, x2max/2, n2x)';   % column vector

%% ========================================================================
%  Part 6: H3 plane (observation plane) parameters
%% ========================================================================

n3x   = 153;                                        % number of sample points on H3 plane
x3max = (n3x - 1) * 3 / scale * 1e-3 / 3;         % half-width of H3 plane (m)
d3x   = x3max / (n3x - 1);                         % sample spacing on H3 plane
ds3   = d3x;                                        % integration element length on H3 plane (for power integration)

% H3-plane x and z coordinates (rotated then translated to absolute coordinate system)
x3_local = -linspace(-x3max/2, x3max/2, n3x);     % local coordinates
x3_abs   = cos_2t * x3_local - x3_shift;           % absolute x coordinates
z3_abs   = -sin_2t * x3_local + z3_shift;          % absolute z coordinates

%% ========================================================================
%  Part 7: Read the reflector surface file and construct H2-plane shape
%% ========================================================================

surf_file = 'divider_surface.txt';
if ~isfile(surf_file)
    error('Reflector surface file not found: %s\nPlease ensure it is in the same directory as this script.', surf_file);
end

% Read node data (35 control points, 1 row)
surf_nodes = dlmread(surf_file);
fprintf('Successfully read reflector surface file: %s\n', surf_file);
fprintf('Number of nodes: %d\n', numel(surf_nodes));

% Back up the loaded surface data
dlmwrite('Divider_surface_backup.txt', surf_nodes, 'delimiter', ' ');

% Interpolate 35 control nodes to n2x=103 sample points
node_idx    = 1 : 3 : 103;              % positions of control nodes within 1~103
interp_idx  = 1 : 103;                  % target interpolation positions
surf_interp = interp1(node_idx, surf_nodes, interp_idx, 'spline');

% Paraboloid baseline shape (focal length: focus)
parab_base = zeros(1, n2x);
for q = -n2x_half : n2x_half
    parab_base(n2x_half + 1 - q) = -((q * d2x)^2) / (4 * focus);
end

% H2-plane total shape = paraboloid baseline + surface perturbation
surf_z2 = parab_base + surf_interp;    % H2-plane z values (local coordinates)

%% ========================================================================
%  Part 8: Compute rotated coordinates of H2 plane (absolute coordinate system)
%% ========================================================================

x2_rot = zeros(1, n2x);    % absolute x coordinate of each H2-plane point
z2_rot = zeros(1, n2x);    % absolute z coordinate after rotation

for i = 1 : n2x
    x2_rot(i) =  cos_t * x2_coords(i) + sin_t * surf_z2(i);
    z2_rot(i) = -sin_t * x2_coords(i) + cos_t * surf_z2(i);
end
z2_abs = z2_shift + z2_rot;    % add z-direction translation to get absolute z coordinates

%% ========================================================================
%  Part 9: Compute the normal vectors of H2-plane points (after rotation)
%  Used for direction cosine (obliquity factor) computation
%% ========================================================================

% Finely interpolate surf_z2 (half step size) for numerical derivative estimation of normals
interp_fine_idx  = 1 : 0.5 : 103;
surf_z2_fine = interp1(interp_idx, surf_z2, interp_fine_idx, 'spline');

% Normal vectors in local coordinate system (before rotation)
nx2_local = zeros(1, n2x);
nz2_local = -ones(1, n2x);             % default normal pointing downward (-z)

for i = 2 : n2x - 1
    dz = surf_z2_fine(2*i-1) - surf_z2_fine(2*i-2);  % z-direction finite difference
    dx = d2x;                                          % x-direction step size
    norm_len = sqrt(dz^2 + dx^2);
    nx2_local(i) =  dz / norm_len;
    nz2_local(i) = -dx / norm_len;
end

% Rotate to absolute coordinate system
nx2_abs = zeros(1, n2x);
nz2_abs = zeros(1, n2x);
for i = 1 : n2x
    nx2_abs(i) =  cos_t * nx2_local(i) + sin_t * nz2_local(i);
    nz2_abs(i) = -sin_t * nx2_local(i) + cos_t * nz2_local(i);
end

%% ========================================================================
%  Part 10: First propagation — H1 plane -> H2 plane
%  2D scalar diffraction integral (Helmholtz-Kirchhoff)
%% ========================================================================

% Horn aperture x coordinates (arranged right-to-left, consistent with original code)
x1_horn = zeros(1, n1x);
for i = 1 : n1x
    x1_horn(i) = horn_aperture/2 - (i-1) * horn_aperture / (n1x-1);
end

% H2-plane magnetic field (computed by superimposing contributions from H1 plane)
H2 = zeros(1, n2x);
cos_factor12 = zeros(n2x, n1x);

for i = 1 : n2x
    h_contrib = zeros(1, n1x);
    for j = 1 : n1x
        r12 = sqrt(z2_abs(i)^2 + (x1_horn(j) - x2_rot(i))^2);
        % Direction cosine (obliquity factor) = cos(angle between observation direction and z-axis)
        cos_factor12(i, j) = z2_abs(i) / r12;
        % 2D Helmholtz integral kernel: Hankel function H0^(2)(kr)
        h_contrib(j) = ds1 * H1(j) * cos_factor12(i,j) ...
                       * besselh(0, 2, beta * r12) * pi / lambda;
    end
    H2(i) = sum(h_contrib);
end

%% ========================================================================
%  Part 11: Second propagation — H2 plane -> H3 plane
%  H2 plane acts as an equivalent surface source radiating to H3
%% ========================================================================

H3 = zeros(1, n3x);
h3_contrib = zeros(1, n2x);

for u3 = 1 : n3x
    for u2 = 1 : n2x
        rxx = x3_abs(u3) - x2_rot(u2);
        rzz = z3_abs(u3) - z2_abs(u2);
        rll = sqrt(rxx^2 + rzz^2);
        % Direction cosine (including normal vector projection)
        cos_factor23 = (rzz * nz2_abs(u2) + rxx * nx2_abs(u2)) / rll;
        h3_contrib(u2) = cos_factor23 * H2(u2) ...
                         * besselh(0, 2, beta * rll) * pi / lambda;
    end
    H3(u3) = sum(h3_contrib) * ds2;
end

%% ========================================================================
%  Part 12: Output results
%% ========================================================================

% --- Output power ---
P_output = 100 * 0.5 * 377 * sum(abs(H3).^2) * ds3 * 2;
fprintf('H3 plane normalized output power P_output = %.6f\n', P_output);

% --- Magnitude (dB, clipped to -30 dB) ---
H3_abs  = abs(H3);
H3_max  = max(H3_abs);
H3_dB   = 20 * log10(H3_abs / H3_max);
H3_dB(H3_dB < -30) = -30;             % dynamic range clipping

% --- Phase (degrees) ---
Phase_offset = 100;                    % phase offset (degrees), consistent with original code
H3_phase = angle(H3 * exp(j_unit * Phase_offset * pi / 180)) * 180 / pi;

% --- Print fitness metrics ---
H3_left  = H3(1 : floor(n3x/2));
H3_right = H3(ceil(n3x/2)+1 : end);
power_left  = sum(abs(H3_left).^2)  * ds3;
power_right = sum(abs(H3_right).^2) * ds3;
balance = abs(power_left - power_right) / (power_left + power_right);
fprintf('Left half power  = %.6f\n', power_left);
fprintf('Right half power = %.6f\n', power_right);
fprintf('Power imbalance  = %.4f%%\n', balance * 100);

% --- Plot ---
figure(1);
subplot(2, 1, 1);
plot(H3_dB, 'b', 'LineWidth', 2);
xlabel('H3 Plane Sample Index');
ylabel('Magnitude (dB)');
title('H3 Plane Magnetic Field Magnitude Distribution');
grid on;
ylim([-30, 0]);

subplot(2, 1, 2);
plot(H3_phase, 'r', 'LineWidth', 2);
xlabel('H3 Plane Sample Index');
ylabel('Phase (deg)');
title('H3 Plane Magnetic Field Phase Distribution');
grid on;

% Overlay reference curves if target files exist
if isfile('Divider_mag.txt')
    ref_mag = dlmread('Divider_mag.txt');
    subplot(2, 1, 1); hold on;
    plot(ref_mag, 'r--', 'LineWidth', 1.5);
    legend('Computed', 'Reference');
end
if isfile('Divider_phase.txt')
    ref_phase = dlmread('Divider_phase.txt');
    subplot(2, 1, 2); hold on;
    plot(ref_phase, 'r--', 'LineWidth', 1.5);
    legend('Computed', 'Reference');
end

% --- Correlation coefficient with target field pattern (if reference file exists) ---
if isfile('Divider_C.txt')
    scale_factor = 0.45;
    H3_ref = scale_factor * dlmread('Divider_C.txt');
    H3_shifted = H3 * exp(j_unit * (Phase_offset - 1) * pi / 180);
    corr_coeff = sum(H3_shifted .* conj(H3_ref)) ...
                 / sqrt(sum(abs(H3_shifted).^2) * sum(abs(H3_ref).^2));
    fprintf('Correlation coefficient with target field k = %.6f\n', abs(corr_coeff));
end

toc;
fprintf('\nComputation complete.\n');
