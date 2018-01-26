% Define parameters/variables
M = sym('M');
M1 = sym('M1');
M2 = sym('M2');
L1 = sym('L1');
L2 = sym('L2');
t = sym('t');
theta1 = sym('theta1');
ttheta1 = sym('theta1(t)');
theta2 = sym('theta2');
ttheta2 = sym('theta2(t)');
dtheta1 = sym('dtheta1');
tdtheta1 = sym('tdtheta1(t)');
dtheta2 = sym('dtheta2');
tdtheta2 = sym('tdtheta2(t)');
q = sym('q');
tq = sym('q(t)');
dq = sym('dq');
tdq = sym('tdq(t)');
g = sym('g');
u = sym('u');


% Kinetic energy
K = sym('0.5*M*dq^2 + 0.5*M1*((dq + L1*dtheta1*cos(theta1))^2 + (L1*dtheta1*sin(theta1))^2) + 0.5*M2*((dq+L1*dtheta1*cos(theta1)+L2*dtheta2*cos(theta1+theta2))^2+(L1*dtheta1*sin(theta1)+L2*dtheta2*sin(theta1+theta2))^2)'); 

% Potential energy
P = sym('M1*g*L1*cos(theta1)+M2*g*(L1*cos(theta1)+L2*cos(theta1+theta2))'); 

% Lagrangian
L(q,dq,theta1,theta2,dtheta1,dtheta2) = K-P;

% Differentiation
dLdq = diff(L,q)
dLddq = diff(L,dq)
dLdtheta1 = diff(L,theta1)
dLdtheta2 = diff(L,theta2)
dLddtheta1 = diff(L,dtheta1)
dLddtheta2 = diff(L,dtheta2)

ddtdLdq = diff(dLdq(tq,tdq,ttheta1,ttheta2,tdtheta1,tdtheta2),t)
ddtdLddq = diff(dLddq(tq,tdq,ttheta1,ttheta2,tdtheta1,tdtheta2),t)
ddtdLdtheta1 = diff(dLdtheta1(tq,tdq,ttheta1,ttheta2,tdtheta1,tdtheta2),t)
ddtdLdtheta2 = diff(dLdtheta1(tq,tdq,ttheta1,ttheta2,tdtheta1,tdtheta2),t)
ddtdLddtheta1 = diff(dLddtheta1(tq,tdq,ttheta1,ttheta2,tdtheta1,tdtheta2),t)
ddtdLddtheta2 = diff(dLddtheta2(tq,tdq,ttheta1,ttheta2,tdtheta1,tdtheta2),t)

E1 = ddtdLddq-dLdq
E2 = ddtdLddtheta1 - dLdtheta1
E3 = ddtdLddtheta2 - dLdtheta2

