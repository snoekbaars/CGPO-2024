System_full;

Rk = 1e-5*[2.5*eye(3) zeros(3); 
           zeros(3) 7.57*eye(3)];

Qk = [[10 0 0;
     0 10 0;
     0 0 1000] zeros(3,9); zeros(9,3) eye(9)];
B1 = eye(12);
%Qk = B1*Q*B1';
Qk = 1e2*eye(4);


%ksys = ss(Ad, [Bd zeros(12); zeros(4) eye(4)], Cd, [Dd zeros(6, 4)]); % Fout
[kalmf, L, P] = kalman(sysd, Qk, Rk);

System;
LQI;
