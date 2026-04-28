%==========================================================================
% generate_reflector_stl.m
%
% Function: Generate the 3D shape of a power divider reflector surface by
%           superimposing external perturbations onto a standard paraboloid,
%           and export it as an STL file for use in EM simulation software
%           (e.g., HFSS).
%
% Workflow:
%   1. Build the x/y coordinate grid of the reflector surface
%   2. Compute the standard paraboloid baseline shape
%   3. Read perturbation data from divider_surface.txt and interpolate
%   4. Superimpose the paraboloid and perturbation to obtain the final shape
%   5. Preview plot and export STL file
%
% Input file:
%   divider_surface.txt  — reflector perturbation node data (1x35, unit: m)
%
% Output file:
%   reflector_surface.stl  — reflector triangular mesh (unit: mm, ASCII format)
%
%==========================================================================

clc;
clear;

%% ========================================================================
%  Part 1: Basic parameter settings
%% ========================================================================

% Reference scaling factor (normalization coefficient for the design baseline frequency)
scale = 3.8;

% Paraboloid focal length (m)
focal_length = 0.11 / scale;

% Number of sample points per segment between control nodes (interpolation density)
pts_per_segment = 3;

%% ========================================================================
%  Part 2: x-axis coordinate settings (reflector width direction)
%% ========================================================================

n_x      = 103;                                          % number of sample points in x direction
x_span   = (n_x - 1) * 3 / scale * 1e-3 / pts_per_segment;  % total width in x direction (m)
dx       = x_span / (n_x - 1);                          % sample spacing in x direction (m)
n_x_half = (n_x - 1) / 2;                               % number of points from center to edge

% x coordinate array (arranged from positive to negative, consistent with EM coordinate system)
x_coords = -linspace(-x_span/2, x_span/2, n_x)';       % column vector (m)

%% ========================================================================
%  Part 3: y-axis coordinate settings (reflector thickness direction,
%           corresponding to waveguide narrow wall)
%% ========================================================================

n_y    = 3;                          % number of sample points in y direction (thin-plate structure, 3 layers suffice)
y_span = 0.254e-3;                   % total thickness in y direction (m), corresponding to standard waveguide narrow wall dimension

% y coordinate array (arranged from positive to negative)
y_coords = -linspace(-y_span/2, y_span/2, n_y)';       % column vector (m)

%% ========================================================================
%  Part 4: Compute standard paraboloid baseline shape
%  Paraboloid equation: z = -x^2 / (4f), where f is the focal length
%% ========================================================================

parab_z = zeros(1, n_x);
for q = -n_x_half : n_x_half
    idx = n_x_half + 1 - q;
    parab_z(idx) = -((q * dx)^2) / (4 * focal_length);
end
% parab_z: z values of the standard paraboloid at each x sample point (m)

%% ========================================================================
%  Part 5: Read perturbation data and interpolate to full resolution
%% ========================================================================

surf_file = 'divider_surface.txt';
if ~isfile(surf_file)
    error('Reflector perturbation file not found: %s\nPlease ensure it is in the same directory as this script.', surf_file);
end

% Read perturbation values of 35 control nodes (m)
surf_ctrl_pts = dlmread(surf_file);
fprintf('Successfully read perturbation file: %s, number of control nodes: %d\n', surf_file, numel(surf_ctrl_pts));

% Interpolate 35 control nodes (one every 3 points) to 103 sample points
ctrl_node_idx = 1 : 3 : 103;    % indices of control nodes within 1~103
full_idx      = 1 : 103;         % target full-resolution indices
surf_perturb  = interp1(ctrl_node_idx, surf_ctrl_pts, full_idx, 'spline');
% surf_perturb: interpolated perturbation values, 1x103 (m)

%% ========================================================================
%  Part 6: Superimpose paraboloid and perturbation to obtain the final
%           reflector shape
%% ========================================================================

% Final reflector z values = standard paraboloid + perturbation
reflector_z_1d = parab_z + surf_perturb;   % 1xn_x, unit: m

% Replicate along y direction to form a 2D surface (n_y x n_x)
reflector_z_2d = repmat(reflector_z_1d, n_y, 1);

%% ========================================================================
%  Part 7: Preview plot
%% ========================================================================

figure(1);

% Subplot 1: final reflector profile (x-z cross-section)
subplot(2, 2, 1);
plot(x_coords * 1e3, reflector_z_1d * 1e3, 'b', 'LineWidth', 1.5);
xlabel('x (mm)');
ylabel('z (mm)');
title('Reflector Profile (with perturbation)');
grid on;

% Subplot 2: 3D shape of the reflector
subplot(2, 2, 2);
surfl(x_coords * 1e3, y_coords * 1e3, reflector_z_2d * 1e3);
xlabel('x (mm)');
ylabel('y (mm)');
zlabel('z (mm)');
title('Reflector 3D Shape');
shading interp;

% Subplot 3: perturbation distribution
subplot(2, 2, 3);
plot(x_coords * 1e3, surf_perturb * 1e6, 'r', 'LineWidth', 1.5);
xlabel('x (mm)');
ylabel('Perturbation (um)');
title('Reflector Perturbation Distribution');
grid on;

% Subplot 4: standard paraboloid profile (for reference comparison)
subplot(2, 2, 4);
plot(x_coords * 1e3, parab_z * 1e3, 'k--', 'LineWidth', 1.5);
xlabel('x (mm)');
ylabel('z (mm)');
title('Standard Paraboloid Baseline (no perturbation)');
grid on;

%% ========================================================================
%  Part 8: Export STL file (unit converted to mm)
%% ========================================================================

stl_filename = 'reflector_surface.stl';

% stlwrite requires coordinates in mm, so multiply by 1000
stlwrite(stl_filename, ...
         x_coords * 1e3, ...
         y_coords * 1e3, ...
         reflector_z_2d * 1e3, ...
         'mode', 'ascii');

fprintf('STL file exported: %s\n', stl_filename);
fprintf('Done.\n');
