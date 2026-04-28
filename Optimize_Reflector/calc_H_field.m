function H3 = calc_H_field(surf_nodes)
%==========================================================================
% calc_H_field  —— Electromagnetic field computation function (called by the optimization loop)
%
% Function: Given the perturbation values of the reflector control nodes,
%           compute the magnetic field distribution on the H3 observation
%           plane of a terahertz power divider using the Physical Optics
%           (PO) method.
%
% Computation flow:
%   H1 (horn aperture) -> free-space propagation -> H2 (reflector) -> free-space propagation -> H3 (observation plane)
%
% Input:
%   surf_nodes  — 1x35 real vector, perturbation values of 35 reflector
%                 control nodes (unit: m)
%
% Output:
%   H3          — 1x153 complex vector, H-field distribution on H3 plane
%
% Note: This function contains no file I/O or plotting. All parameters are
%       defined internally and are consistent with compute_H_from_surface.m.
%==========================================================================

%% ========================================================================
%  Part 1: Basic physical parameters
%% ========================================================================

freq_GHz = 380;
lambda   = 300 / freq_GHz * 1e-3;   % free-space wavelength (m)
beta     = 2 * pi / lambda;          % phase constant (rad/m)
j_unit   = sqrt(-1);

scale = 3.8;   % reference scaling factor

%% ========================================================================
%  Part 2: Geometric parameters
%% ========================================================================

theta  = 45 * pi / 180;             % reflector tilt angle (rad)
cos_t  = cos(theta);
sin_t  = sin(theta);
cos_2t = cos(2 * theta);
sin_2t = sin(2 * theta);

focus    = 0.11 / scale;            % paraboloid focal length (m)
z2_shift = 66 / 3.6 * 1e-3;        % H2-plane z-axis translation (m)
z_HFSS   = 10 / scale * 1e-3;      % HFSS reference plane offset (m)
z3_shift = (z2_shift * sin_2t * tan(2*theta) - z_HFSS) * cos_2t;
x3_shift = 50 / scale * 1e-3;      % H3-plane x-axis translation (m)

%% ========================================================================
%  Part 3: Horn aperture (H1 plane) parameters
%% ========================================================================

horn_length   = 22 / scale * 1e-3;
wg_width      = 0.508e-3;
horn_aperture = 20 / scale * 1e-3;
R_horn        = horn_length / (1 - (wg_width / horn_aperture));

n1x = 15;
ds1 = horn_aperture / (n1x - 1);
x1s = linspace(-horn_aperture/2, horn_aperture/2, n1x);

%% ========================================================================
%  Part 4: H2 plane (reflector surface) parameters
%% ========================================================================

n2x      = 103;
pts_per  = 3;
x2max    = (n2x - 1) * 3 / scale * 1e-3 / pts_per;
d2x      = x2max / (n2x - 1);
ds2      = d2x;
n2x_half = (n2x - 1) / 2;

x2s = -linspace(-x2max/2, x2max/2, n2x)';   % column vector (m)

%% ========================================================================
%  Part 5: H3 plane (observation plane) parameters
%% ========================================================================

n3x   = 153;
x3max = (n3x - 1) * 3 / scale * 1e-3 / pts_per;

x3_local = -linspace(-x3max/2, x3max/2, n3x);
x3_abs   = cos_2t * x3_local - x3_shift;    % global x coordinates of H3 plane
z3_abs   = -sin_2t * x3_local + z3_shift;   % global z coordinates of H3 plane

%% ========================================================================
%  Part 6: Compute H-field distribution at H1 plane (TE10-mode horn aperture field)
%% ========================================================================

E0 = 377;
H1 = zeros(1, n1x);
for idx = 1 : n1x
    x = x1s(idx);
    % Cosine amplitude distribution x spherical-wave phase (equivalent phase center at horn vertex)
    H1(idx) = -E0 / 377 * cos(pi * x / horn_aperture) ...
              * exp(-j_unit * beta * x^2 / (2 * R_horn));
end

%% ========================================================================
%  Part 7: Build H2-plane shape (paraboloid baseline + interpolated perturbation)
%% ========================================================================

% Paraboloid baseline shape
parab_z = zeros(1, n2x);
for q = -n2x_half : n2x_half
    parab_z(n2x_half + 1 - q) = -((q * d2x)^2) / (4 * focus);
end

% Interpolate 35 control nodes (every 3 points) to 103 sample points
ctrl_idx    = 1 : 3 : 103;
full_idx    = 1 : 103;
surf_interp = interp1(ctrl_idx, surf_nodes, full_idx, 'spline');

% Total z values of H2 plane (local coordinate system)
surf_z2 = parab_z + surf_interp;

%% ========================================================================
%  Part 8: Rotate H2-plane coordinates to global coordinate system
%% ========================================================================

x2_rot = zeros(1, n2x);
z2_rot = zeros(1, n2x);
for i = 1 : n2x
    x2_rot(i) =  cos_t * x2s(i) + sin_t * surf_z2(i);
    z2_rot(i) = -sin_t * x2s(i) + cos_t * surf_z2(i);
end
z2_abs = z2_shift + z2_rot;

%% ========================================================================
%  Part 9: Compute H2-plane normal vectors (after rotation, for direction cosines)
%% ========================================================================

% Finely interpolate surf_z2 (step 0.5) for numerical derivative
fine_idx     = 1 : 0.5 : 103;
surf_z2_fine = interp1(full_idx, surf_z2, fine_idx, 'spline');

nx2_local = zeros(1, n2x);
nz2_local = -ones(1, n2x);

for i = 2 : n2x - 1
    dz = surf_z2_fine(2*i-1) - surf_z2_fine(2*i-2);
    len = sqrt(dz^2 + d2x^2);
    nx2_local(i) =  dz  / len;
    nz2_local(i) = -d2x / len;
end

% Rotate to global coordinate system
nx2_abs = zeros(1, n2x);
nz2_abs = zeros(1, n2x);
for i = 1 : n2x
    nx2_abs(i) =  cos_t * nx2_local(i) + sin_t * nz2_local(i);
    nz2_abs(i) = -sin_t * nx2_local(i) + cos_t * nz2_local(i);
end

%% ========================================================================
%  Part 10: H1 -> H2 propagation (2D Helmholtz-Kirchhoff diffraction integral)
%  H2(i) = sum_j [ H1(j) * (z2/r12) * H0^(2)(beta*r12) * pi/lambda * ds1 ]
%% ========================================================================

H2 = zeros(1, n2x);
for i = 1 : n2x
    h_sum = 0;
    for j = 1 : n1x
        r12 = sqrt(z2_abs(i)^2 + (x1s(j) - x2_rot(i))^2);
        cos_factor12 = z2_abs(i) / r12;
        h_sum = h_sum + ds1 * H1(j) * cos_factor12 * besselh(0, 2, beta * r12) * pi / lambda;
    end
    H2(i) = h_sum;
end

%% ========================================================================
%  Part 11: H2 -> H3 propagation (equivalent surface source radiation with normal vector projection)
%  H3(u) = sum_i [ H2(i) * (rz*nz2 + rx*nx2)/r23 * H0^(2)(beta*r23) * pi/lambda * ds2 ]
%% ========================================================================

H3 = zeros(1, n3x);
for u = 1 : n3x
    h_sum = 0;
    for i = 1 : n2x
        rx  = x3_abs(u) - x2_rot(i);
        rz  = z3_abs(u) - z2_abs(i);
        r23 = sqrt(rx^2 + rz^2);
        cos_factor23 = (rz * nz2_abs(i) + rx * nx2_abs(i)) / r23;
        h_sum = h_sum + ds2 * H2(i) * cos_factor23 * besselh(0, 2, beta * r23) * pi / lambda;
    end
    H3(u) = h_sum;
end

end
