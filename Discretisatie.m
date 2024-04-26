Ts = 0.05;

A = zeros(12);
B = zeros(12,4);
C = zeros(6, 12);
D = zeros(6, 4);

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
