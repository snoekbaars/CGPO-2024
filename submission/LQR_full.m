%% LQR state-feedback controller of form (Nx, Nu, K)
clear all; close all; clc;

% Loading the discretized linearized system
System_full;

% The system is stored in the variable sysd
sysd

Q = diag([10 10 1000 1 1 1 1 1 1 1 1 1]);

R = eye(4);

[K, S, P] = lqr(sysd, Q, R);

% The number of inputs p is not equal to the number of outputs
% so we need to use the pseudo-inverse function 

N = pinv([Ad - eye(12) Bd; Cd  Dd])*[zeros(12, 3); eye(3); zeros(3)];
Nx = N(1:12,:);
Nu = N(13:end, :);

