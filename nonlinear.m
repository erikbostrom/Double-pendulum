%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverted double pendulum
% nonlinear model around equilibrium in origo.
% Erik Bostrom, erikbos@kth.se 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nonlinear()
% Parameters
g=9.81;
l1 = 2;
l2 = 1;
m = 100;
m1 = 10;
m2 = 10;

% Initial conditions
epsilon = 0.01;
x0 = [-0.5 0 0 0 0 0]';
xe0 = [0.5 0 0 0 0 0]';
xhat0 = x0;
%xhat0 = x0 + epsilon;
%xhat0 = [1 0 0 0 0 0]';

K =   [21
      -11618
       34166
          69
       -1663
        8664]';
    
L = [24.0000    0         0
    -11.9896   11.9294    0.0602
    11.9500  -35.9106   23.9606
    128.9810   -0.4905   -0.4905
   -133.5520  199.0939  -65.5419
     72.9876 -278.7294  205.7418 ];



%% Solution
t0=0; tend=10; tend2=2;

% Ode solver
[T, X] =ode45(@xdot,[t0 tend],x0,[],g,l1,l2,m,m1,m2,K);
[T2, X2] =ode45(@xdot2,[t0 tend2],xe0,[],g,l1,l2,m,m1,m2,L);
[T3, X3] =ode45(@xdot3,[t0 tend],[x0 xhat0],[],g,l1,l2,m,m1,m2,K,L);


%% Plots

 figure(1); plot(T,X(:,1),'k-'); xlabel('time'); ylabel('q');   
 figure(2); plot(T,X(:,2),'k-'); xlabel('time'); ylabel('\theta_1');
 figure(3); plot(T,X(:,3),'k-'); xlabel('time'); ylabel('\theta_2');
 figure(4); plot(T2,X2(:,1),'k-'); xlabel('time'); ylabel('e_1');
 figure(5); plot(T2,X2(:,2),'k-'); xlabel('time'); ylabel('e_2');
 figure(6); plot(T2,X2(:,3),'k-'); xlabel('time'); ylabel('e_3');
 figure(7); plot(T3,X3(:,1),'k-'); xlabel('time'); ylabel('q');
 figure(8); plot(T3,X3(:,2),'k-'); xlabel('time'); ylabel('\theta_1');
 figure(9); plot(T3,X3(:,3),'k-'); xlabel('time'); ylabel('\theta_2');
return


%% Systems of first order differential equations
function dxdt= xdot(t,x,g,l1,l2,m,m1,m2,K)
a11 = m+m1+m2;
a12 = l1*(m1+m2)*cos(x(2));
a13 = l2*m2*cos(x(2)+x(3));
a21 = l1*(m1+m2)*cos(x(2));
a22 = l1^2*(m1+m2);
a23 = l1*l2*m2*cos(x(3));
a31 = l2*m2*cos(x(2)+x(3));
a32 = l1*l2*m2*cos(x(3));
a33 = l2^2*m2;
%K = 1.0e+04 * [0.0021   -1.1618    3.4166    0.0069   -0.1663    0.8664];
u = -K*x;
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

Kappa = cond(M)
return

% Error dynamics edot = f(t,e,u=0)-L*C*e
function dxdt= xdot2(t,x,g,l1,l2,m,m1,m2,L)
a11 = m+m1+m2;
a12 = l1*(m1+m2)*cos(x(2));
a13 = l2*m2*cos(x(2)+x(3));
a21 = l1*(m1+m2)*cos(x(2));
a22 = l1^2*(m1+m2);
a23 = l1*l2*m2*cos(x(3));
a31 = l2*m2*cos(x(2)+x(3));
a32 = l1*l2*m2*cos(x(3));
a33 = l2^2*m2;

%L=[24.0000   -0.0000    0.0000
%  -11.9896   11.9294    0.0602
%   11.9500  -35.9106   23.9606
%  128.9810   -0.4905   -0.4905
% -133.5520  199.0939  -65.5419
%   72.9876 -278.7294  205.7418];

 C = [1   0    0  0  0  0;
      1   l1   0  0  0  0;
      1 l1+l2 l2  0  0  0];
 
M = [a11 a12 a13;
     a21 a22 a23;
     a31 a32 a33];

f=[l1*(m1+m2)*x(5)^2*sin(x(2)) + l2*m2*x(5)*x(6)*sin(x(2)+x(3))+l2*m2*x(6)^2*sin(x(2)+x(3));% + u;
   g*l1*(m1+m2)*sin(x(2))+g*l2*m2*sin(x(2)+x(3))-l2*m2*x(4)*x(6)*sin(x(2)+x(3)) + l1*l2*m2*x(6)^2*sin(x(3));
   g*l2*m2*sin(x(2)+x(3))+l2*m2*x(4)*x(5)*sin(x(2)+x(3))];

dxdt(1,1) = x(4);
dxdt(2,1) = x(5);
dxdt(3,1) = x(6);
dxdt(4:6,1) = M\f;

dxdt(:,1) = dxdt(:,1) - L*C*[x(1);x(2);x(3);x(4);x(5);x(6)];
return
 
% Combined state feedback and state observer
function dxdt= xdot3(t,x,g,l1,l2,m,m1,m2,K,L)

a11 = m+m1+m2;
a12 = l1*(m1+m2)*cos(x(2));
a13 = l2*m2*cos(x(2)+x(3));
a21 = l1*(m1+m2)*cos(x(2));
a22 = l1^2*(m1+m2);
a23 = l1*l2*m2*cos(x(3));
a31 = l2*m2*cos(x(2)+x(3));
a32 = l1*l2*m2*cos(x(3));
a33 = l2^2*m2;

%K = 1.0e+04 * [0.0021   -1.1618    3.4166    0.0069   -0.1663    0.8664];

%L=[24.0000   -0.0000    0.0000
%  -11.9896   11.9294    0.0602
%   11.9500  -35.9106   23.9606
%  128.9810   -0.4905   -0.4905
% -133.5520  199.0939  -65.5419
%   72.9876 -278.7294  205.7418];
 
M = [a11 a12 a13;
     a21 a22 a23;
     a31 a32 a33];

u=-K*x(7:12);

f=[l1*(m1+m2)*x(5)^2*sin(x(2)) + l2*m2*x(5)*x(6)*sin(x(2)+x(3))+l2*m2*x(6)^2*sin(x(2)+x(3)) + u;
   g*l1*(m1+m2)*sin(x(2))+g*l2*m2*sin(x(2)+x(3))-l2*m2*x(4)*x(6)*sin(x(2)+x(3)) + l1*l2*m2*x(6)^2*sin(x(3));
   g*l2*m2*sin(x(2)+x(3))+l2*m2*x(4)*x(5)*sin(x(2)+x(3))];

f2=[l1*(m1+m2)*x(11)^2*sin(x(8)) + l2*m2*x(11)*x(12)*sin(x(8)+x(9))+l2*m2*x(12)^2*sin(x(8)+x(9)) + u;
    g*l1*(m1+m2)*sin(x(8))+g*l2*m2*sin(x(8)+x(9))-l2*m2*x(10)*x(12)*sin(x(8)+x(9)) + l1*l2*m2*x(12)^2*sin(x(9));
    g*l2*m2*sin(x(8)+x(9))+l2*m2*x(10)*x(11)*sin(x(8)+x(9))];

y=zeros(3,1);
y(1)=x(1);
y(2)=x(1)+l1*sin(x(2));
y(3)=x(1)+l1*sin(x(2))+l2*sin(x(2)+x(3));
y1=zeros(3,1);
y1(1)=x(7);
y1(2)=x(7)+l1*sin(x(8));
y1(3)=x(7)+l1*sin(x(8))+l2*sin(x(8)+x(9));

dxdt(1,1) = x(4);
dxdt(2,1) = x(5);
dxdt(3,1) = x(6);
dxdt(4:6,1) = M\f;
dxdt(7,1)=x(10);
dxdt(8,1)=x(11);
dxdt(9,1)=x(12);

dxdt(10:12,1)=M\f2;
dxdt(7:12,1) = dxdt(7:12,1) + L*(y-y1);
return


