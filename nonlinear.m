%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverted double pendulum
% nonlinear model around equilibrium in origo.
% Erik Bostrom, erikbos@kth.se 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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