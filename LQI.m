%% LQR controller with integral action

% Augmenting the system with extra states integrating the output error:

Aa = [eye(3) Cd; zeros(12, 3) Ad];

Ba = [Dd;Bd]; 

Ca = [zeros(3, 3) Cd]; 

Da = Dd; 

sysa = ss(Aa, Ba, Ca, Da, Ts);

% Verify controllability of the system:

Co_a = ctrb(sysa);
rank(Co_a) % Has to be 15, and it is

xyz_weights = 1e2*[1 0 0;
                0 1 0;
                0 0 1;];

Qa = [1e2*eye(3) zeros(3, 12); 
      zeros(3, 3) xyz_weights zeros(3, 9);
      zeros(9, 6) 1e2*eye(9)];


Ra = eye(4);

[Ka, Sa, Pa] = lqi(sysd, Qa, Ra);
K0_a = Ka(:, 1:12);
K1_a = Ka(:, 13:end);