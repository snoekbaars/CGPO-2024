Ts = 0.05;

A = zeros(12);
B = zeros(12,4);
C = zeros(6, 12);
D = zeros(6, 4);

A_d = (eye(12)-A*Ts/2)\(eye(12)+A*Ts/2);
B_d = (eye(12)-A*Ts/2)\B*Ts;
C_d = C/(eye(12)-A*Ts/2);
D_d = D+C/inv(eye(12)-A*Ts/2)*B*Ts/2;

sys = ss(A_d, B_d, C_d, D_d, Ts);
tzero(sys)
pole(sys)
Co = ctrb(sys);
rank(Co); % Moet 12 zijn
Ob = obsv(sys);
rank(Ob)
