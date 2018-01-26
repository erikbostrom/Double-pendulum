g=9.81;
l1 = 2;
l2 = 1;
m = 100;
m1 = 10;
m2 = 10;

x = [0,pi/2,pi/2,0,0,0];

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

MInv = inv(M);

f=[l1*(m1+m2)*x(4)^2 - l2*m2*x(5)*x(6)*sin(x(2)+x(3))+l2*m2*sin(x(2)+x(3));
   g*l1*(m1+m2)*sin(x(2))+g*l2*m2*sin(x(2)+x(3))-l2*m2*x(4)*x(6)*sin(x(2)+x(3)) + l1*l2*m2*x(6)*sin(x(3));
   g*l2*m2*sin(x(2)+x(3))+l2*m2*x(4)*x(5)*sin(x(2)+x(3))];