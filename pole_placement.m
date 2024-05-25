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

[K, S, P] = lqr(sysd, Q, R);

% The number of inputs p is not equal to the number of outputs
% so we need to use the pseudo-inverse function 

N = pinv([Ad - eye(12) Bd; Cd  Dd])*[zeros(12, 3); eye(3); zeros(3)];
Nx = N(1:12,:);
Nu = N(13:end, :);

eig(sysd.A - sysd.B*K)

%% Step 1) defining the desired location of the poles (first without payload)
% There are twelve poles to place

% Smaller (more negative) nondominant poles lead to slower response times but also smaller
% control actions (and larger coefficients in the gain matrix)

% Method: apply pole placement method in continuous time and discretize the
% poles using z = e^(s*Ts)

% Without payload: 

dr = 0.9;  % damping ratio
t_set = 4; % settling time
wn = 4.6/(dr*t_set); % natural frequency
alpha = -dr*wn % real part dominant poles
beta = wn*sqrt(1-dr^2); % imag part dominant poles
% Dominant poles
p_d = [alpha + beta*1i, alpha - beta*1i];
% Non dominant poles
p_nd = [-2.1, -2.11, -2.12, -2.13, -2.14, -2.15, -2.16, -2.17, -2.18, -2.19];

% The continuous time poles are mapped to discrete time poles using e^(s*Ts)
poles_c = exp(Ts*[p_d, p_nd])';

% Using the matlab function place to find the corresponding gain matrix L
Kp = place(sysd.A, sysd.B, poles_c);
Acl = sysd.A - sysd.B*Kp;

% Checking the so-obtained closed loop system and its poles 
syscl = ss(Acl, sysd.B, sysd.C, sysd.D, Ts);
pole(syscl);

% With the above poles selection, the quadcopter goes through all of the
% checkpoints in time.

%% With payload
% When adding a payload, the above pole placement leads to a steady state
% error in the Z-coordinate (the weight of the quadcopter assumed when
% deriving the first-principle laws of the model is no longer valid)

% Needed to increase the non-dominant poles in absolute value, to be able
% to compensate for the added load in order to get rid of the steady state
% error in the Z-dimension. This leads to control actions which exceed the
% limits specified in the assignment:

dr = 0.9;  % damping ratio
t_set = 4; % settling time
wn = 4.6/(dr*t_set); % natural frequency

alpha = -dr*wn % real part dominant poles
beta = wn*sqrt(1-dr^2); % imag part dominant poles

p_d = [alpha + beta*1i, alpha - beta*1i];
p_nd = 2.8*[-2.1, -2.11, -2.12, -2.13, -2.14, -2.15, -2.16, -2.17, -2.18, -2.19];

poles_c = exp(Ts*[p_d, p_nd])';

Kp = place(sysd.A, sysd.B, poles_c);

%% Step 2) Design of the estimator

% We proceed in the same way as for the controller, but this time using the
% rule of thumb that the estimator should have poles 2 to 5 times faster
% than that of the controller

dr = 0.9;  % high damping ratio 
t_set = 1; % settling time - higher than that for the controller

dr = 0.9;  % high damping ratio 
t_set = 1.2; % settling time - higher than that for the controller

dr = 0.9;  % high damping ratio 
t_set = 0.1; % settling time - higher than that for the controller
wn = 4.6/(dr*t_set); % natural frequency
alpha = -dr*wn % real part dominant poles
beta = wn*sqrt(1-dr^2); % imag part dominant poles
% Dominant poles
p_d = [alpha + beta*1i, alpha - beta*1i];
% Non dominant poles
p_nd = [-2.1, -2.11, -2.12, -2.13, -2.14, -2.15, -2.16, -2.17, -2.18, -2.19];
p_nd = [-20, -20.1, -20.11, -20.12, -20.123, -20.2, -20.21, -20.23, -20.3, -20.32];

%p_nd = 0.8*[-20, -20.1, -20.11, -20.12, -20.123, -20.2, -20.21, -20.23, -20.3, -20.32];

p_nd = 1.5*[-20, -21, -22, -23, -24, -25, -26, -27, -28, -29];

% Discrete time poles
poles_c = exp(Ts*[p_d, p_nd])';

% Applying the duality relationship between controller and estimator design

L = place(sysd.A', sysd.C', poles_c)';

Ae = sysd.A-L*sysd.C;
Be = [sysd.B L];
Ce = eye(12);
De = zeros(12, 10);

%%
% We proceed in the same way as for the controller, but this time using the
% rule of thumb that the estimator should have poles 2 to 5 times faster
% than that of the controller

dr = 0.9;  % high damping ratio 
t_set = 1; % settling time - higher than that for the controller

dr = 0.9;  % high damping ratio 
t_set = 1.2; % settling time - higher than that for the controller

dr = 0.9;  % high damping ratio 
t_set = 0.1; % settling time - higher than that for the controller
wn = 4.6/(dr*t_set); % natural frequency
alpha = -dr*wn % real part dominant poles
beta = wn*sqrt(1-dr^2); % imag part dominant poles
% Dominant poles
p_d = [alpha + beta*1i, alpha - beta*1i];
% Non dominant poles
p_nd = [-2.1, -2.11, -2.12, -2.13, -2.14, -2.15, -2.16, -2.17, -2.18, -2.19];
p_nd = [-20, -20.1, -20.11, -20.12, -20.123, -20.2, -20.21, -20.23, -20.3, -20.32];

%p_nd = 0.8*[-20, -20.1, -20.11, -20.12, -20.123, -20.2, -20.21, -20.23, -20.3, -20.32];

p_nd = 1.5*[-20, -21, -22, -23, -24, -25, -26, -27, -28, -29];

% Discrete time poles
poles_c = exp(Ts*[p_d, p_nd])'; 

% poles_c = [0.9456 + 0.0000i
%   -0.0296 + 0.2492i
%   -0.0296 - 0.2492i
%    0.9456 + 0.0000i
%    0.9512 + 0.0000i
%    0.0025 + 0.0000i
%    0.0025 + 0.0000i
%    0.0075 + 0.0000i
%    0.0075 + 0.0000i
%    0.0075 + 0.0000i
%    0.9512 + 0.0000i
%    0.9512 + 0.0000i];

% Applying the duality relationship between controller and estimator design

L = place(sysd.A', sysd.C', poles_c)';

Ae = sysd.A-L*sysd.C;
Be = [sysd.B L];
Ce = eye(12);
De = zeros(12, 10);

%% 
dr = 0.85;  % damping ratio
t_set = 1; % settling time
wn = 4.6/(dr*t_set); % natural frequency
alpha = -dr*wn % real part dominant poles
beta = wn*sqrt(1-dr^2); % imag part dominant poles
% Dominant poles
p_d = [alpha + beta*1i, alpha - beta*1i];
% Non dominant poles
p_nd = [-2.1, -2.11, -2.12, -2.13, -2.14, -2.15, -2.16, -2.17, -2.18, -2.19];

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

%% 

poles_continuous  = [-2.3+1.113940841i,-2.3-1.113940841i,-2.6,-2.6,-2.6,-2.6,-2.8,-2.8,-2.8,-2.8,-3,-3];

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