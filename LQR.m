clear; close all; clc;
addpath("..")

% CONSTANTS DECLARATION
m = 0.5;        % quadcopter mass
L = 0.25;       % quadcopter radius
k = 3e-6;       % propeller lift coefficient
b = 1e-7;       % propeller drag coefficient
g = 9.81;       % acceleration due to gravity
k_d = 0.25;     % air friction coefficient
Ixx = 5e-3;     % inertia about x_b axis
Iyy = 5e-3;     % inertia about y_b axis
Izz = 1e-2;     % inertia about z_b axis
c_m = 1e4;      % motor constant

% 4.1) LINEARIZATION
% Finding the value of the equilibrium point
v1_2 = (1/4) * ((m*g)/(k*c_m));

% input voltages v_1, ... , v_4
v1 = sqrt(v1_2);
v2 = v1;
v3 = v1;
v4 = v1;

% Linearizing the nonlinear model about the equilibrium point

%state_vars = [x y z v_x, v_y, v_z, phi, theta, psi, w_x, w_y, w_z]

% A is 12 by 12
A = [
    0 0 0 1 0 0 0 0 0 0 0 0                 %1
    0 0 0 0 1 0 0 0 0 0 0 0;                %2
    0 0 0 0 0 1 0 0 0 0 0 0;                %3
    0 0 0 (-k_d/m) 0 0 0 g 0 0 0 0;         %4
    0 0 0 0 (-k_d/m) 0 -g 0 0 0 0 0;        %5
    0 0 0 0 0 (-k_d/m) 0 0 0 0 0 0;         %6
    0 0 0 0 0 0 0 0 0 1 0 0;                %7
    0 0 0 0 0 0 0 0 0 0 1 0;                %8
    0 0 0 0 0 0 0 0 0 0 0 1;                %9
    0 0 0 0 0 0 0 0 0 0 0 0;                %10
    0 0 0 0 0 0 0 0 0 0 0 0;                %11
    0 0 0 0 0 0 0 0 0 0 0 0;                %12
    ];

% Input vars = [v1_2, v2_2, v3_2, v4_2]

% B is 12 by 4
B = [
    0 0 0 0;                                            %1  
    0 0 0 0;                                            %2
    0 0 0 0;                                            %3
    0 0 0 0;                                            %4
    0 0 0 0;                                            %5
    (k*c_m)/m (k*c_m)/m (k*c_m)/m (k*c_m)/m;            %6
    0 0 0 0;                                            %7
    0 0 0 0;                                            %8
    0 0 0 0;                                            %9
    (L*k*c_m)/Ixx 0 -(L*k*c_m/Ixx) 0;                   %10
    0 (L*k*c_m)/Iyy 0 -(L*k*c_m)/Iyy;                   %11
    (b*c_m)/Izz -(b*c_m)/Izz (b*c_m)/Izz -(b*c_m)/Izz;  %12
    ];

% C is 6 by 12
C = [
    1 0 0 0 0 0 0 0 0 0 0 0;
    0 1 0 0 0 0 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 1 0 0 0 0 0;
    0 0 0 0 0 0 0 1 0 0 0 0;
    0 0 0 0 0 0 0 0 1 0 0 0;
    ];
% C is 3 by 12
C = [eye(3,3) zeros(3,9)];
D = zeros(6, 4);
D = zeros(3, 4);

% 4.2) DISCRETIZATION - Bilinear
Ts = 0.05;

% A = zeros(12);
% B = zeros(12,4);
% C = zeros(6, 12);
% D = zeros(6, 4);

Ad = (eye(12)-A*Ts/2)\(eye(12)+A*Ts/2);
Bd = (eye(12)-A*Ts/2)\B*Ts;
Cd = C/(eye(12)-A*Ts/2);
Dd = D+C/inv(eye(12)-A*Ts/2)*B*Ts/2;



sys = ss(Ad, Bd, Cd, Dd, Ts);

tzero(sys)
pole(sys)
Co = ctrb(sys);
rank(Co) % Moet 12 zijn
[Vc, Ec] = eig(Ad');
% PBH test
for i = 1:12
    q = Vc(:,i);
    if norm(Bd'*q) < 1e-10
        i
        Ec(i,i)
    end
end
Ob = obsv(sys);
rank(Ob)
[Vo, Eo] = eig(Ad);
for i = 1:12
    p = Vo(:,i);
    if norm(Cd*p) < 1e-10
        i
        Eo(i,i)
    end
end

Q = [[10 0 0;
     0 10 0;
     0 0 10] zeros(3,9); zeros(9,3) eye(9)];
R = 10*eye(4);
[K,S,P] = lqr(sys,Q,R);

load("../references_05.mat")

Ns = pinv([A B; C D])*[zeros(12,3); eye(3)];
Nx = Ns(1:12,:);
Nu = Ns(13:16, :);

