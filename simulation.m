%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Double pendulum 
%
% Time integration of the nonlinear (unstabilized) double 
% pendulum equations.
%
% Dynamical system derived using the Euler--Lagrange equations.
%
% Erik Boström, erikbos@kth.se
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;

% Parameters
g=9.81;
l1 = 1/3;
l2 = 2/3;
m1 = 10;
m2 = 1;
t0=0; 
T=10;
Nt = 1000;

% Initial conditions
x10 = pi+0.1;
x20 = pi+0.1;
x30 = 0;
x40 = 0;


%% Ode solver
It = linspace(t0,T,Nt);
[T, X] =ode45(@xdot_free,It,[x10 x20 x30 x40],[],g,l1,l2,m1,m2);

%% Plots
fig=0;

% Plot solution
figure(fig+1)
subplot(121)
plot(T,X(:,3),'k-');
xlabel('time','Interpreter','latex'); ylabel('$\theta_1$','Interpreter','latex');
    
subplot(122)
plot(T,X(:,4),'k-');
xlabel('time','Interpreter','latex'); ylabel('$\theta_2$','Interpreter','latex');
saveas(gcf,'imgs/theta.png');
close all;

% Movie
x1=l1*sin(X(:,1));
y1=-l1*cos(X(:,1));
x2=l1*sin(X(:,1))+l2*sin(X(:,2));
y2=-l1*cos(X(:,1))-l2*cos(X(:,2));

filename = 'imgs/simulation.gif';
h=figure(fig+1);
str = sprintf('Double pendulum: $l_1=%.1f$, $l_2=%.1f$, $m_1=%.1f$, $m_2=%.1f$',l1,l2,m1,m2);
for j = 1:Nt
    plot([0 x1(j)],[0 y1(j)],'k-o',[x1(j) x2(j)],[y1(j) y2(j)],'k-o');
    axis([-1.2*(l1+l2) 1.2*(l1+l2) -1.2*(l1+l2) 1.2*(l1+l2)]);    
    title(str, 'Interpreter', 'latex');
    xlabel('$x$', 'Interpreter', 'latex'); ylabel('$y$', 'Interpreter', 'latex');
    drawnow;
    frame = getframe(h);
    im    = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    % Write to gif movie
    if j == 1 
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0); 
    else 
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0); 
    end 
end

% Plot trajectories
figure(fig+1)
plot(x1,y1,'r-',x2,y2,'k-',[0],[0],'k*');
axis([-1.2*(l1+l2) 1.2*(l1+l2) -1.2*(l1+l2) 1.2*(l1+l2)]);
axis square;
title('Trajectories of the pendulums','Interpreter', 'latex');
xlabel('$x$','Interpreter', 'latex'); ylabel('$y$','Interpreter', 'latex');
legend('#1','#2','Center axis');
saveas(gcf,'imgs/trajectory.png');




