close all;
clear variables;
clc;

addpath images\ dependencies\;

%% Parameters of the vehicle

c1 = 1500; % Damping coefficients (Ns/m)
c2 = c1;

k1 = 21000; % Stiffness coefficients (N/m)
k2 = k1;

kt1 = 150000; % Tire stiffness coefficients (N/m)
kt2 = kt1;

m1 = 40; % Unsprung masses (kg)
m2 = m1;

m = 400; % Sprung mass (kg)

Iyy = 600; % Inertia (kg/m^2)
ar = 1.45; % Distance from CG to rear-axle (m)
af = 0.8; % Distance from CG to front-axle (m)
L_ = ar+af; % Distance between front-axle and rear-axle (m)

V = 20/3.6; % Speed (m/s)
dt = L_/V; % Delay between front and rear

%% State space representation of the system

A1 = [-c1/m-af^2*c1/Iyy c1/m+af^2*c1/Iyy -c2/m+af*ar*c2/Iyy c2/m-af*ar*c2/Iyy;
    c1/m1 -c1/m1 0 0;
    -c1/m+af*ar*c1/Iyy c1/m-af*ar*c1/Iyy -c2/m-ar^2*c2/Iyy c2/m-ar^2*c2/Iyy;
    0 0 c2/m2 -c2/m2];

A2 = [-k1/m-af^2*k1/Iyy k1/m+af^2*k1/Iyy -k2/m+af*ar*k2/Iyy k2/m-af*ar*k2/Iyy;
    k1/m1 -(k1+kt1)/m1 0 0;
    -k1/m+af*ar*k1/Iyy k1/m-af*ar*k1/Iyy -k2/m-ar^2*k2/Iyy k2/m+ar^2*k2/Iyy;
    0 0 k2/m2 -(k2+kt2)/m2];

A3 = eye(4);
A4 = zeros(4);

A = [A1 A2;
    A3 A4]; 

clear A1 A2 A3 A4;

B1 = [1/m+af^2/Iyy 1/m-af*ar/Iyy;
    -1/m1 0;
    1/m-af*ar/Iyy 1/m+ar^2/Iyy;
    0 -1/m2];

B2 = zeros(4,2);

B = [B1;
    B2];

clear B1 B2;

K1 = [0 0;
    kt1/m1 0;
    0 0;
    0 kt2/m2];

K2 = zeros(4,2);

K = [K1;
    K2];

clear K1 K2;

C = [0 0 0 0 1 0 0 0;
     0 0 0 0 0 0 1 0];

D = 0;

% P = idss(A, B, C, D, K);
P = ss(A, B, C, D); % Creating the state-space entity

%% H-infinity loop-shaping controller design

om_rad = logspace(-4, 4, 501);
bw_rad = 30; % Target BW

figure(1)
sigma(P, om_rad)
legend('P')
grid on;

%% Get the system to the desired bandwidth

[singVals, ~] = sigma(P, bw_rad);
maxSingVals = max(singVals);

Wo1 = eye(2)/maxSingVals;

figure(1)
clf 
sigma(P, 'k--', Wo1*P, 'b-', om_rad)
legend('P', 'Wo1*P')
grid on;

%% Add in a PI weight to improve disturbance rejection and avoid 
%% Steady state error

Om_pi = 1*bw_rad; % Apply the weight at the target frequency

W_pi = tf([1 Om_pi], [1 0]);

Wo2 = Wo1*mdiag(W_pi, 10*W_pi);

figure(1)
clf
sigma(P, 'r:', Wo1*P, 'k--', Wo2*P, 'b-', om_rad)
legend('P', 'Wo1*P', 'Wo2*P')
grid on;

%% Add in a HF filter to improve noise rejection

Om_hf = 1*bw_rad; % Apply the weight at the target frequency

W_hf = tf(Om_hf, [1 Om_hf]);

Wo3 = Wo2*mdiag(W_hf, W_hf);

figure(1)
clf
sigma(P, 'r:', Wo1*P, 'k--', Wo2*P, 'b--', Wo3*P, 'b-', om_rad)
legend('P', 'Wo1*P', 'Wo2*P', 'Wo3*P')
grid on;

Wo = Wo3; % Output weight which is a combination of Proporttional, Pi and HF
Wi = eye(2); % Input weight untouched

Pw = Wo*P*Wi; % Weighted plant

Cinf_minus = ncfsyn(Pw); % Hinf synthesis
Cinf = -Cinf_minus;

figure(1)
clf
sigma(P, 'r:', Pw, 'k--', Pw*Cinf, 'b-', om_rad)
legend('P', 'Pw', 'Pw*Cinf')
grid on;


C = Wi*Cinf*Wo; % Our controller

%% Characteristics of the system when the Controller is applied

L = P*C;
T = feedback(L, eye(2));
S = eye(2) - T;
CS = feedback(C, P);
SP = feedback(P, C);

figure(2)
clf
subplot(121)
rho(Pw, Cinf_minus)
subplot(222)
bodemag(S, 'k-', T, 'b-', CS, 'k:')
grid on;
legend('S', 'T', 'CS', 'Location','southwest')
subplot(224)
step(T, 'b', 10)
grid on;