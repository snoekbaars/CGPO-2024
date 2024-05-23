%% LQR controller with integral action

clear all; close all; clc;

% Loading the discretized linearized system
System_full;

% The system is stored in the variable sysd
sysd

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