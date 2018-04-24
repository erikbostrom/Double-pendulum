%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamical system for the stabilized inverted double pendulum problem
%
% t  : time
% x  = (x1,...x6) : state variable vector
%   where:
%   x1 = theta1
%   x2 = theta2
%   x3 = dx1/dt
%   x4 = dx2/dt
%
% l1 : length of inner pendulum
% l2 : length of outer pendulum
% m1 : mass of weight at middle joint
% m2 : mass of outer weight
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dxdt = xdot_free(t,x,g,l1,l2,m1,m2)

dxdt(1,1) = x(3);
dxdt(2,1) = x(4);
dxdt(3,1) = ( -m2*cos(x(1)+x(2))*l1*sin(x(1)-x(2))*x(3)^2 ...
              +m2*cos(x(1)-x(2))*g*sin(x(2)) ...
              -m2*l2*sin(x(1)-x(2))*x(4)^2 ...
              -(m1+m2)*g*sin(x(1)) )...
              /(l1*(m1+m2)-m2*cos(x(1)-x(2))^2);
dxdt(4,1) = ( (m1+m2)*(l1*sin(x(1)-x(2))*x(3)^2) ...
    + (sin(x(1)-x(2))*cos(x(1)-x(2))*m2*l2*x(4)^2)/(m1+m2) ...
    + cos(x(1)+x(2))*g*sin(x(1))-g*sin(x(2)) )...
    /(l2*(m1+m2*sin(x(1)+x(2))^2));
return
