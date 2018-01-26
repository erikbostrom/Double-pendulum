% %--------------------------------------------------------------
% Dynamical Systems
%
% Inverted double pendulum
% Solution to the non-linear problem
% By Erik Boström
%--------------------------------------------------------------

function nonlinear2
close all;

% Parameters
g=9.81;
l1 = 1;
l2 = 1;
m = 100;
m1 = 10;
m2 = 10;

% Initial conditions
x10 = 0;
x20 = pi/4;
x30 = pi/4;
x40 = 0;
x50 = 0;
x60 = 0;
u0 = 0;


%% Solution
t0=0; tend=10;

% Ode solver
[T, X] =ode45(@xdot,[t0:0.02:tend],[x10 x20 x30 x40 x50 x60],[],u0,g,l1,l2,m,m1,m2);

%% Plots
plot(T,X(:,5),'k-');
title('Non-linear solution');
xlabel('time'); ylabel('\theta_1');
    
figure(2);
plot(T,X(:,6),'k-');
title('Non-linear solution');
xlabel('time'); ylabel('\theta_2');
return


%% System of first order differential equations
function dxdt= xdot(t,x,u,g,l1,l2,m,m1,m2)
a11 = m+m1+m2;
a12 = l1*(m1+m2)*cos(x(2));
a13 = l2*m2*cos(x(2)+x(3));
a21 = l1*(m1+m2)*cos(x(2));
a22 = l1^2*(m1+m2);
a23 = l1*l2*m2*cos(x(3));
a31 = l2*m2*cos(x(2)+x(3));
a32 = l1*l2*m2*cos(x(3));
a33 = l2^2*m2;

M = [a11 a12 a13;
     a21 a22 a23;
     a31 a32 a33];

f=[l1*(m1+m2)*x(5)^2*sin(x(2)) + l2*m2*x(5)*x(6)*sin(x(2)+x(3))+l2*m2*x(6)^2*sin(x(2)+x(3)) + u;
   g*l1*(m1+m2)*sin(x(2))+g*l2*m2*sin(x(2)+x(3))-l2*m2*x(4)*x(6)*sin(x(2)+x(3)) + l1*l2*m2*x(6)^2*sin(x(3));
   g*l2*m2*sin(x(2)+x(3))+l2*m2*x(4)*x(5)*sin(x(2)+x(3))];

dxdt(1,1) = x(4);
dxdt(2,1) = x(5);
dxdt(3,1) = x(6);
dxdt(4:6,1) = M\f;
return

