%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverted double pendulum
% Linearized problem around the upper equilibruim
% (u=q=theta1=theta2=0)
%
% Erik Bostrom, erikbos@kth.se 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
g=9.81;
L1 = 2;
L2 = 1;
M = 100;
M1 = 10;
M2 = 10;

a0 = -g*(L1*(M1 + M2) + L2*M2)/(L1*M);
a1 = -g*L2*M2/(L1*M);
a2 = g*(L1*M1*(M+M1+M2)+L2*M2*(M+M1))/(L1^2*M*M1);
a3 = g*M2*(-L1*M+L2*(M+M1))/(L1^2*M*M1);
a4 = -g*M2/(L1*M1);
a5 = g*(L1*(M1+M2)-L2*M2)/(L1*L2*M1);

%% State space model
A = [0  0  0  1  0  0;
     0  0  0  0  1  0;
     0  0  0  0  0  1;
     0 a0 a1  0  0  0;
     0 a2 a3  0  0  0;
     0 a4 a5  0  0  0];

B = [0 0 0 1/M -1/(L1*M) 0]';

C = [1   0    0  0  0  0;
     1   L1   0  0  0  0;
     1 L1+L2 L2  0  0  0];

C2 = [1   0   0  0  0  0;
     0   1   0  0  0  0;
     0   0   1  0  0  0];
 
D = zeros(3,1);

% State space form
sys = ss(A,B,C,[]);

% Eigenvalues of A
l = eig(A);

% Pole Placement
p = [-1 -2 -1+1i -1-1i -2+1i -2-1i];
K = place(A,B,p);
%K = acker(A,B,p2);

sys2 = ss(A-B*K,B,C2,D);

% State observer
pl = 8*p;
L=place(A',C',pl).';
est = estim(sys,L);

% Controlled system
Ae=[A       -B*K     ;
    L*C   A-L*C-B*K];
Be = [B; B];
cosys = ss(Ae,Be,[],[]);

% Controllability matrix
R = ctrb(sys);

% Observability matrix
Ob = obsv(sys);

%% Plots
...