%% Linear Quadratic Gaussian Control 
% Using the LQR with integral action and 
% a Kalman filter for the state estimation

clear all; close all; clc;

% Loading the discretized linearized system
System_full;

% The system is stored in the variable sysd
sysd

%% LQR with Integral Action

% Matlab has a built-in function lqi for computing the LQR
% with integral action. This function assumes that the dimensions
% of the output and the reference are the same. We therefore truncate
% the matrices C and D of the original system
C_lqi = Cd(1:3, :);
D_lqi = Dd(1:3, :);

% New system: 

sys_lqi = ss(Ad, Bd, C_lqi, D_lqi, Ts);

% Augmenting the system with extra states integrating the output

Aa = [eye(3) C_lqi; zeros(12, 3) Ad];

Ba = [D_lqi;Bd]; 

Ca = [zeros(3, 3) C_lqi]; 

Da = D_lqi; 

% Creating augmented system to test for the controllability
sysa = ss(Aa, Ba, Ca, Da, Ts);
Co_a = ctrb(sysa);
rank(Co_a) % Rank = 15, as required

% Computing the lqi

% Setting matrices Q and R for augmented system
Qa = 1e2*diag([1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]);

Ra = eye(4);

% Computing the optimal gain
[Ka, Sa, Pa] = lqi(sys_lqi, Qa, Ra);
K_s = Ka(:, 1:12);
K_i = Ka(:, 13:end);

%% Designing the Kalman filter

% We are using the kalman function from matlab to design the kalman filter

% To use this function, we need to take into account the matrix B1 provided
% and create a new system, with a modified B matrix. This will be the
% system with noise.

Bk = [sysd.B eye(12)];
Dk = [sysd.D zeros(6, 12)];

sysk = ss(sysd.A, Bk, sysd.C, Dk, Ts);

% Estimate of the variance of the process noise
var_p = 1e-2;

% Process noise covariance matrix with weighting
Qk = var_p*diag([1 1 1e-6 1 1 20 1 1 1 1 1 1]);

% Measurement noise
Rk = diag([2.5e-5 2.5e-5 2.5e-5 7.57e-5 7.57e-5 7.57e-5]);

[kalmf, L, P] = kalman(sysk, Qk, Rk);

%% Estimator gain

% The estimation error (without noise) is governed by the system matrix
syse_matrix = (sysd.A - L*sysd.C);

eig(syse_matrix)



