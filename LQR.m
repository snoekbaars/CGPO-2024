Q = [[10 0 0;
     0 10 0;
     0 0 1000] zeros(3,9); zeros(9,3) eye(9)];
R = 0.5*eye(4);

[K,S,P] = lqr(sysd,Q,R);

Ns = pinv([Ad-eye(12) Bd; Cd Dd])*[zeros(12,3); eye(3)];
Nx = Ns(1:12,:);
Nu = Ns(13:16, :);