%--------------------------------------------------------------
% Dynamical Systems
%
% Inverted double pendulum
% Linearized problem around the equilibruim
% u=q=theta1=theta2=0
% By Erik Boström
%--------------------------------------------------------------

close all;
fprintf('Linearized problem of the inverted double pendulum around\nthe upper equilibrum point.\n\n');

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

B = [0 0 0 1/M1 -1/(L1*M) 0]';

C = [1   0    0  0  0  0;
     1   L1   0  0  0  0;
     1 L1+L2 L2  0  0  0];

D = zeros(3,1);

% State space form
sys = ss(A,B,C,D);
 
% Eigenvalues of A
l = eig(A);

% Pole Placement
p = [-1 -1+1i -1-1i -2 -2+1i -2-1i];
K = place(A,B,p);
sys2 = ss(A-B*K,B,C,D);

ll = eig(A-B*K); % Check that eigenvalues are indeed p

% Controllability matrix
R = ctrb(sys);
fprintf('The controllability matrix: R=\n'); disp(R);
fprintf('rank R = %d\n\n',rank(R));

% Observability matrix
Ob = obsv(sys);
fprintf('The observability matrix: O=\n\n'); disp(Ob);

% Transient responses
x0=[-5 0 0 0 0 0]';
step(A,B,C,D,[],'k');
xlabel( 't' );
ylabel( '\theta');
laprint(1,'step1')

figure(2)
step(A-B*K,B,C,D,[],'k');

figure(11)
initial(sys,x0,'k');
figure(12)
initial(sys2,x0,'k');
S = stepinfo(sys2);


%% Plots

% Bode plot
figure(3);
bode(sys,'k');

figure(4)
bode(sys2,'k');

% Eigenvalues
figure(5);
plot(l,[0 0 0 0 0 0],'ko');
axis square; grid on;
xlabel('Re'); ylabel('Im');
title('Eigenvalues');
axis([-5 5 -5 5]);
figure(6);
pzmap(sys2,'k');
title('Eigenvalues controlled system')
axis([-3 3 -3 3]);