%% LQR state-feedback controller
clear all; close all; clc;

% Loading the discretized linearized system
System;
% The system is stored in the variable sysd
sysd

Q = diag([10 10 1000 1 1 1 1 1 1 1 1 1]);

R = 1*eye(4);

[K,S,P] = lqr(sysd,Q,R);

Ns = pinv([Ad-eye(12) Bd; Cd Dd])*[zeros(12,3); eye(3)];
Nx = Ns(1:12,:);
Nu = Ns(13:16, :);