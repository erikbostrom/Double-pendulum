function test2
close all;

l1 = 2;
l2 = 1;

%% Solution
t0=0; tend=8;
x10 = 0;
x20 = pi/2;
x30 = 0;
x40 = 0;
x50 = 0;
x60 = 0;
u0=0;

[T, X] =ode45(@dx456,[t0:0.01:tend],[x10 x20 x30 x40 x50 x60],u0);

%% Plots
plot(T,X(:,5),'k-');
title('Non-linear solution');
xlabel('time'); ylabel('\theta_1');

figure(2);
plot(T,X(:,6),'k-');
title('Non-linear solution');
xlabel('time'); ylabel('\theta_2');

x1=L1*sin(X(:,1));
y1=-L1*cos(X(:,1));
x2=L1*sin(X(:,1))+L2*sin(X(:,2));
y2=-L1*cos(X(:,1))-L2*cos(X(:,2));

plot(x1,y1);
figure(2);
plot(x2,y2);
return


function [xp] = dx456(t,x,u)
x1=x(1);
x2=x(2);
x3=x(3);
x4=x(4);
x5=x(5);
x6=x(6);

c2 = cos(x2);
s2 = sin(x2);
c3 = cos(x3);
s3 = sin(x3);
c23 = cos(x2+x3);
s23 = sin(x2+x3);
V = [ 120    40*c2 10*c23;
      40*c2  80    20*c3;
      10*c23 20*c3 10];
W = [10*(x5+x6)*x6*s23 + 40*x5^2*s2 + u;
 (100 - 10*x4*x6)*s23 + 400*s2 + 20*x6^2*s3;
 (100 + 10*x4*x5)*s23 ];
xp = inv(V)*W
return