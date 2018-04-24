%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamical system for the stabilized inverted double pendulum problem
%
% t  : time
% x  = (x1,...x6) : state variable vector
%   where:
%   x1 = u
%   x2 = theta1
%   x3 = theta2
%   x4 = dx1/dt
%   x5 = dx2/dt
%   x6 = dx3/dt
%
% u  : force on cart
% l1 : length of inner pendulum
% l2 : length of outer pendulum
% m  : mass of cart
% m1 : mass of weight at middle joint
% m2 : mass of outer weight
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dxdt = xdot(t,x,u,g,l1,l2,m,m1,m2)
a11 = m+m1+m2;
a12 = l1*(m1+m2)*cos(x(2));
a13 = l2*m2*cos(x(2)+x(3));
a21 = l1*(m1+m2)*cos(x(2));
a22 = l1^2*(m1+m2);
a23 = l1*l2*m2*cos(x(3));
a31 = l2*m2*cos(x(2)+x(3));
a32 = l1*l2*m2*cos(x(3));
a33 = l2^2*m2;

M = [ a11 a12 a13 ; 
      a21 a22 a23 ;
      a31 a32 a33 ];

f1 = l1*(m1+m2)*x(5)^2*sin(x(2)) + l2*m2*x(5)*x(6)*sin(x(2)+x(3)) + ...
     l2*m2*x(6)^2*sin(x(2)+x(3)) + u;
f2 = g*l1*(m1+m2)*sin(x(2)) + g*l2*m2*sin(x(2)+x(3)) - ...
     l2*m2*x(4)*x(6)*sin(x(2)+x(3)) + l1*l2*m2*x(6)^2*sin(x(3));
f3 = g*l2*m2*sin(x(2)+x(3)) + l2*m2*x(4)*x(5)*sin(x(2)+x(3));

f = [ f1 ; 
      f2 ; 
      f3 ];

dxdt(1,1) = x(4);
dxdt(2,1) = x(5);
dxdt(3,1) = x(6);
dxdt(4:6,1) = M\f;
return
