%% Normoptimal_allpassv2.m
% Normoptimale ILC with a all pass 
% setup
clear all
close 
s = tf('s');
sys = tf(1,[1.6 1]);
sys2 = tf(1,[1.6 1]);
% sys = (1-2*s)/((s+2)*(s+5));
% sys = (1-s)/(s+2)/(s+3);
%sys = (1-s)/(s-2)/(s-3);
T=.01;
sysd = c2d(sys,T);
sys2d = c2d(sys,T);
[A,b,c,d] = ssdata(sysd);
%% Calculate the first column of G and Toeplitz
    g(1) = c*b;
    M = eye(1);
    for ii = 2:100
       M = M*A;
       g(ii) = c*M*b;
    end
% Generate the  full Toeplitz Matrix
    Gvoll = toeplitz(g);
% Extract the lower trinabular part
    G = tril(Gvoll);
%% Set the weighting matrices and Learn factor L
Wdelta = 0.1*eye(100);
We = eye(100);
% Learn factor
L = inv(Wdelta + G'*We*G) * G' * We;
x0 = [0;0];
% Reference Input and time is 101 long
r = [3*ones(1,20) 1*ones(1,20) 4 *ones(1,20) 2*ones(1,20)  1*ones(1,21) ]';
t = [0:.01:1]';
% Initialize the plant input
uold = zeros(101,1);
u = zeros(101,1);
% Initialize the error
eold = zeros(101,1);
e = zeros(101,1);
% run the iteration loop
imax = 5000; %Number of iteration cycles
for ii = 1:imax    %loop over the iterations
    u(1:100) = uold(1:100) + L * eold(2:101);
    [y,t] = lsim(sysd,u,t,x0,'zoh');
    sum = 0;
    for k = 1:100
        e(k+1) = r(k+1) - y(k+1);
        sum = sum + e(k+1)^2;
    end
    e2(ii) = sum;
    il(ii) = ii;
    eold = e;
    uold = u;
end
% Plot the results of the error and the respsonse after imax iterations
subplot(321),semilogy(il,e2)
%xlabel('$j$','interpreter','latex'),
ylabel('$||e(k)||^2_2$','interpreter','latex')
strwdelta=num2str(Wdelta(1,1)); strwe=num2str(We(1,1));
title(['{\bf Optimization Based:} ${\bf W}_{\Delta u}$  = ' strwdelta ' ${\bf W}_e$ = ' strwe],'interpreter','latex')
axis([0 5000 0.01 1e2])
subplot(322),plot(t/T,y,'g',t/T,r,'--r'),xlabel('$k$','interpreter','latex')
ylabel('- $y(k)$, - - - $r(k)$','interpreter','latex')
axis([0 100 0 5])
%% Compare with the gradient method - Owen07a
% Initialize the plant input
uold = zeros(101,1);
u = zeros(101,1);
% Initialize the error
eold = zeros(101,1);
e = zeros(101,1);
% for a safe model beta<2
beta = 20;
% beta =0.8;
%
imax = 5000; %Number of iteration cycles
for ii = 1:imax    %loop over the iterations
    u(1:100) = uold(1:100) +  beta * G' * eold(2:101);
    u(101) = 0;  % not needed relative degree r=1
    [y,t] = lsim(sysd,u,t,x0,'zoh');
    sum = 0;
    for k = 2:101
        e(k) = r(k) - y(k);
        sum = sum + e(k)^2;
    end
    sum21(ii) = sum;
    zyklus(ii) = ii;
    eold = e;
    uold = u;
end
% Plot the results of the error and the respsonse after imax iterations
subplot(323),semilogy(zyklus,sum21)
%xlabel('$j$','interpreter','latex'),
ylabel('$||e(k)||^2_2$','interpreter','latex')
strbeta = num2str(beta);
title(['Gradient Based: \beta = ' strbeta])
axis([0 5000 1e-5 1e2])
subplot(324),plot(t/T,y,'g',t/T,r,'--r'),xlabel('$k$','interpreter','latex')
ylabel('- $y(k)$, - - - $r(k)$','interpreter','latex')
axis([0 100 0 5])
%% Compare with Inversion based
% Initialize the plant input
uold = zeros(101,1);
u = zeros(101,1);
% Initialize the error
eold = zeros(101,1);
e = zeros(101,1);
% for a safe model beta<2
beta = .25;
beta = 1.8;
%
clear y zyklus sum21
imax = 200; %Number of iteration cycles
for ii = 1:imax    %loop over the iterations
    u(1:100) = uold(1:100) +  beta * inv(G) * eold(2:101);
    u(101) = 0;   % not needed relative degree r=1
    [y,t] = lsim(sysd,u,t,x0,'zoh');
    sum = 0;
    for k = 2:101
        e(k) = r(k) - y(k);
        sum = sum + e(k)^2;
    end
    sum21(ii) = sum;
    zyklus(ii) = ii;
    eold = e;
    uold = u;
end
% Plot the results of the error and the respsonse after imax iterations
subplot(325),semilogy(zyklus,sum21),hold
xlabel('$j$','interpreter','latex'),ylabel('$||e(k)||^2_2$','interpreter','latex')
strbeta = num2str(beta);
title(['Inversion Based: \beta = ' strbeta])
axis([0 200 1e-30 1e2])
subplot(326),plot(t/T,y,'g',t/T,r,'--r'),xlabel('$k$','interpreter','latex')
ylabel('- $y(k)$, - - - $r(k)$','interpreter','latex')
axis([0 100 0 5])
%print -depsc2 normoptimal_allpass_st1.eps