System_full;

Qk = 1e-6*eye(4);
Rk = 1e-5*[2.5*eye(3) zeros(3); 
           zeros(3) 7.57*eye(3)];

Q = [[10 0 0;
     0 10 0;
     0 0 1000] zeros(3,9); zeros(9,3) eye(9)];
B1 = [zeros(4, 8), eye(4)];
Qk = B1*Q*B1';

[kalmf, L, P] = kalman(sysd, Qk, Rk, zeros(4,6));

System;
LQR;
