%% Compensator design via Pole Placement method
% Clearing workspace 
clear all; close all; clc;

% Loading the discretized linearized system
System_full;

% The system is stored in the variable sysd
sysd

%% LQR based feedback

Q = diag([10 10 1000 1 1 1 1 1 1 1 1 1]);

R = eye(4);

% The number of inputs p is not equal to the number of outputs
% so we need to use the pseudo-inverse function 

N = pinv([Ad - eye(12) Bd; Cd  Dd])*[zeros(12, 3); eye(3); zeros(3)];
Nx = N(1:12,:);
Nu = N(13:end, :);

%% Step 1) defining the desired location of the poles (first without payload)
% There are twelve poles to place

% Smaller (more negative) nondominant poles lead to slower response times but also smaller
% control actions (and larger coefficients in the gain matrix)

% Method: apply pole placement method in continuous time and discretize the
% poles using z = e^(s*Ts)

dr = 0.85;  % damping ratio
t_set = 3; % settling time
wn = 4.6/(dr*t_set); % natural frequency
alpha = -dr*wn; % real part dominant poles
beta = wn*sqrt(1-dr^2); % imag part dominant poles
% Dominant poles
p_d = [alpha + beta*1i, alpha - beta*1i];
% Non dominant poles
p_nd = [-2, -2, -2.2, -2.2, -2.4, -2.4, -2.6, -2.6, -2.7, -2.7];

poles_continuous = [p_d, p_nd];

% The continuous time poles are mapped to discrete time poles using e^(s*Ts)
poles_c = exp(Ts*poles_continuous)';

% Using the matlab function place to find the corresponding gain matrix L
Kp = place(sysd.A, sysd.B, poles_c);

poles_estimator = 5*poles_continuous;

poles_e = exp(Ts*poles_estimator)';

L = place(sysd.A', sysd.C', poles_e)';

Ae = sysd.A-L*sysd.C;
Be = [sysd.B-L*sysd.D L];
Ce = eye(12);
De = zeros(12, 10);