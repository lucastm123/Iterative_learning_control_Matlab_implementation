clc
clear all
close all

% Assignment #4
% With uncertainty


%% system

clear all
close 
s = tf('s');
sys = tf(1,[1 1]);


T = 0.01;
sysd = c2d(sys,T);

[A,b,c,d] = ssdata(sysd);



%%  four extreme case models and plot Nyquist

tfun(1) = tf(1.5,[0.2 1]);
tfun(2) = tf(1.5,[3.0 1]);
tfun(3) = tf(0.5,[0.2 1]);
tfun(4) = tf(0.5,[3.0 1]);
% convert to discrete time using ZOH
tfun_d(1) = c2d( tfun(1) , T);
tfun_d(2) = c2d( tfun(2) , T);
tfun_d(3) = c2d( tfun(3) , T);
tfun_d(4) = c2d( tfun(4) , T);

% figure(1) 
% nyquist(tfun_d(1)/sysd,'b',tfun_d(2)/sysd,'r',tfun_d(3)/sysd,'g',tfun_d(4)/sysd,'k'),hold

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
%% Set the weighting matrices and Learn factor L - 1
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
    [y,t] = lsim(tfun_d(1),u,t,x0,'zoh');
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
figure(1)
% Plot the results of the error and the respsonse after imax iterations
subplot(421),semilogy(il,e2)
%xlabel('$j$','interpreter','latex'),
ylabel('$||e(k)||^2_2$','interpreter','latex')
strwdelta=num2str(Wdelta(1,1)); strwe=num2str(We(1,1));
title(['{\bf Optimization Based:} ${\bf W}_{\Delta u}$  = ' strwdelta ' ${\bf W}_e$ = ' strwe],'interpreter','latex')
axis([0 5000 0.01 1e2])
subplot(422),plot(t/T,y,'g',t/T,r,'--r'),xlabel('$k$','interpreter','latex')
ylabel('- $y(k)$, - - - $r(k)$','interpreter','latex')
axis([0 100 0 5])
title('{\bf $G(s) = \frac{1.5}{(0.2 s + 1)} $}','interpreter','latex')

%% Set the weighting matrices and Learn factor L - 2
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
    [y,t] = lsim(tfun_d(2),u,t,x0,'zoh');
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
figure(1)
% Plot the results of the error and the respsonse after imax iterations
subplot(423),semilogy(il,e2)
%xlabel('$j$','interpreter','latex'),
ylabel('$||e(k)||^2_2$','interpreter','latex')
strwdelta=num2str(Wdelta(1,1)); strwe=num2str(We(1,1));
title(['{\bf Optimization Based:} ${\bf W}_{\Delta u}$  = ' strwdelta ' ${\bf W}_e$ = ' strwe],'interpreter','latex')
axis([0 5000 0.01 1e2])
subplot(424),plot(t/T,y,'g',t/T,r,'--r'),xlabel('$k$','interpreter','latex')
ylabel('- $y(k)$, - - - $r(k)$','interpreter','latex')
axis([0 100 0 5])
title('{\bf $G(s) = \frac{1.5}{(3 s + 1)} $}','interpreter','latex')

%% Set the weighting matrices and Learn factor L - 2
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
    [y,t] = lsim(tfun_d(3),u,t,x0,'zoh');
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
figure(1)
% Plot the results of the error and the respsonse after imax iterations
subplot(425),semilogy(il,e2)
%xlabel('$j$','interpreter','latex'),
ylabel('$||e(k)||^2_2$','interpreter','latex')
strwdelta=num2str(Wdelta(1,1)); strwe=num2str(We(1,1));
title(['{\bf Optimization Based:} ${\bf W}_{\Delta u}$  = ' strwdelta ' ${\bf W}_e$ = ' strwe],'interpreter','latex')
axis([0 5000 0.01 1e2])
subplot(426),plot(t/T,y,'g',t/T,r,'--r'),xlabel('$k$','interpreter','latex')
ylabel('- $y(k)$, - - - $r(k)$','interpreter','latex')
axis([0 100 0 5])
title('{\bf $G(s) = \frac{0.5}{(0.2 s + 1)} $}','interpreter','latex')

%% Set the weighting matrices and Learn factor L - 4
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
    [y,t] = lsim(tfun_d(4),u,t,x0,'zoh');
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
figure(1)
% Plot the results of the error and the respsonse after imax iterations
subplot(427),semilogy(il,e2)
%xlabel('$j$','interpreter','latex'),
ylabel('$||e(k)||^2_2$','interpreter','latex')
strwdelta=num2str(Wdelta(1,1)); strwe=num2str(We(1,1));
title(['{\bf Optimization Based:} ${\bf W}_{\Delta u}$  = ' strwdelta ' ${\bf W}_e$ = ' strwe],'interpreter','latex')
axis([0 5000 0.01 1e2])
subplot(428),plot(t/T,y,'g',t/T,r,'--r'),xlabel('$k$','interpreter','latex')
ylabel('- $y(k)$, - - - $r(k)$','interpreter','latex')
axis([0 100 0 5])
title('{\bf $G(s) = \frac{0.5}{(3 s + 1)} $}','interpreter','latex')

saveas(figure (1),'Q2_Optimization_ILC','epsc')

%%
%% Compare with the gradient method - Owen07a - 1
% Initialize the plant input
uold = zeros(101,1);
u = zeros(101,1);
% Initialize the error
eold = zeros(101,1);
e = zeros(101,1);
% for a safe model beta<2
beta = 1.8;
% beta =0.8;
%
imax = 5000; %Number of iteration cycles
for ii = 1:imax    %loop over the iterations
    u(1:100) = uold(1:100) +  beta * G' * eold(2:101);
    u(101) = 0;  % not needed relative degree r=1
    [y,t] = lsim(tfun_d(1),u,t,x0,'zoh');
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
figure(2)
% Plot the results of the error and the respsonse after imax iterations
subplot(421),semilogy(zyklus,sum21)
%xlabel('$j$','interpreter','latex'),
ylabel('$||e(k)||^2_2$','interpreter','latex')
strbeta = num2str(beta);
title(['Gradient Based: \beta = ' strbeta])
axis([0 5000 1e-5 1e2])
subplot(422),plot(t/T,y,'g',t/T,r,'--r'),xlabel('$k$','interpreter','latex')
ylabel('- $y(k)$, - - - $r(k)$','interpreter','latex')
axis([0 100 0 5])
title('{\bf $G(s) = \frac{1.5}{(0.2 s + 1)} $}','interpreter','latex')

%% Compare with the gradient method - Owen07a - 2
% Initialize the plant input
uold = zeros(101,1);
u = zeros(101,1);
% Initialize the error
eold = zeros(101,1);
e = zeros(101,1);
% for a safe model beta<2
beta = 1.8;
% beta =0.8;
%
imax = 5000; %Number of iteration cycles
for ii = 1:imax    %loop over the iterations
    u(1:100) = uold(1:100) +  beta * G' * eold(2:101);
    u(101) = 0;  % not needed relative degree r=1
    [y,t] = lsim(tfun_d(2),u,t,x0,'zoh');
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
figure(2)
% Plot the results of the error and the respsonse after imax iterations
subplot(423),semilogy(zyklus,sum21)
%xlabel('$j$','interpreter','latex'),
ylabel('$||e(k)||^2_2$','interpreter','latex')
strbeta = num2str(beta);
title(['Gradient Based: \beta = ' strbeta])
axis([0 5000 1e-5 1e2])
subplot(424),plot(t/T,y,'g',t/T,r,'--r'),xlabel('$k$','interpreter','latex')
ylabel('- $y(k)$, - - - $r(k)$','interpreter','latex')
axis([0 100 0 5])
title('{\bf $G(s) = \frac{1.5}{(3 s + 1)} $}','interpreter','latex')

%% Compare with the gradient method - Owen07a - 3
% Initialize the plant input
uold = zeros(101,1);
u = zeros(101,1);
% Initialize the error
eold = zeros(101,1);
e = zeros(101,1);
% for a safe model beta<2
beta = 1.8;
% beta =0.8;
%
imax = 5000; %Number of iteration cycles
for ii = 1:imax    %loop over the iterations
    u(1:100) = uold(1:100) +  beta * G' * eold(2:101);
    u(101) = 0;  % not needed relative degree r=1
    [y,t] = lsim(tfun_d(3),u,t,x0,'zoh');
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
figure(2)
% Plot the results of the error and the respsonse after imax iterations
subplot(425),semilogy(zyklus,sum21)
%xlabel('$j$','interpreter','latex'),
ylabel('$||e(k)||^2_2$','interpreter','latex')
strbeta = num2str(beta);
title(['Gradient Based: \beta = ' strbeta])
axis([0 5000 1e-5 1e2])
subplot(426),plot(t/T,y,'g',t/T,r,'--r'),xlabel('$k$','interpreter','latex')
ylabel('- $y(k)$, - - - $r(k)$','interpreter','latex')
axis([0 100 0 5])
title('{\bf $G(s) = \frac{0.5}{(0.2 s + 1)} $}','interpreter','latex')

%% Compare with the gradient method - Owen07a - 4
% Initialize the plant input
uold = zeros(101,1);
u = zeros(101,1);
% Initialize the error
eold = zeros(101,1);
e = zeros(101,1);
% for a safe model beta<2
beta = 1.8;
% beta =0.8;
%
imax = 5000; %Number of iteration cycles
for ii = 1:imax    %loop over the iterations
    u(1:100) = uold(1:100) +  beta * G' * eold(2:101);
    u(101) = 0;  % not needed relative degree r=1
    [y,t] = lsim(tfun_d(4),u,t,x0,'zoh');
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
figure(2)
% Plot the results of the error and the respsonse after imax iterations
subplot(427),semilogy(zyklus,sum21)
%xlabel('$j$','interpreter','latex'),
ylabel('$||e(k)||^2_2$','interpreter','latex')
strbeta = num2str(beta);
title(['Gradient Based: \beta = ' strbeta])
axis([0 5000 1e-5 1e2])
subplot(428),plot(t/T,y,'g',t/T,r,'--r'),xlabel('$k$','interpreter','latex')
ylabel('- $y(k)$, - - - $r(k)$','interpreter','latex')
axis([0 100 0 5])
title('{\bf $G(s) = \frac{0.5}{(3 s + 1)} $}','interpreter','latex')

saveas(figure (2),'Q2_Gradient_ILC','epsc')

%%
%% Inversion based - 1
% Initialize the plant input
uold = zeros(101,1);
u = zeros(101,1);
% Initialize the error
eold = zeros(101,1);
e = zeros(101,1);
% for a safe model beta<2
beta = 1/3.7;
% beta = 0.8;
%
clear y zyklus sum21
imax = 2000; %Number of iteration cycles
for ii = 1:imax    %loop over the iterations
    u(1:100) = uold(1:100) +  beta * inv(G) * eold(2:101);
    u(101) = 0;   % not needed relative degree r=1
    [y,t] = lsim(tfun_d(1),u,t,x0,'zoh');
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
figure(3)
% Plot the results of the error and the respsonse after imax iterations
subplot(421),semilogy(zyklus,sum21),hold
xlabel('$j$','interpreter','latex'),ylabel('$||e(k)||^2_2$','interpreter','latex')
strbeta = num2str(beta);
title(['Inversion Based: \beta = ' strbeta])
axis([0 200 1e-30 1e2])
subplot(422),plot(t/T,y,'g',t/T,r,'--r'),xlabel('$k$','interpreter','latex')
ylabel('- $y(k)$, - - - $r(k)$','interpreter','latex')
axis([0 100 0 5])
title('{\bf $G(s) = \frac{1.5}{(0.2 s + 1)} $}','interpreter','latex')

%% Inversion based - 2
% Initialize the plant input
uold = zeros(101,1);
u = zeros(101,1);
% Initialize the error
eold = zeros(101,1);
e = zeros(101,1);
% for a safe model beta<2
beta = 1/3.7;
% beta = 0.8;
%
clear y zyklus sum21
imax = 2000; %Number of iteration cycles
for ii = 1:imax    %loop over the iterations
    u(1:100) = uold(1:100) +  beta * inv(G) * eold(2:101);
    u(101) = 0;   % not needed relative degree r=1
    [y,t] = lsim(tfun_d(2),u,t,x0,'zoh');
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
figure(3)
% Plot the results of the error and the respsonse after imax iterations
subplot(423),semilogy(zyklus,sum21),hold
xlabel('$j$','interpreter','latex'),ylabel('$||e(k)||^2_2$','interpreter','latex')
strbeta = num2str(beta);
title(['Inversion Based: \beta = ' strbeta])
axis([0 200 1e-30 1e2])
subplot(424),plot(t/T,y,'g',t/T,r,'--r'),xlabel('$k$','interpreter','latex')
ylabel('- $y(k)$, - - - $r(k)$','interpreter','latex')
axis([0 100 0 5])
title('{\bf $G(s) = \frac{1.5}{(3 s + 1)} $}','interpreter','latex')

%% Inversion based - 3
% Initialize the plant input
uold = zeros(101,1);
u = zeros(101,1);
% Initialize the error
eold = zeros(101,1);
e = zeros(101,1);
% for a safe model beta<2
beta = 1/3.7;
% beta = 0.8;
%
clear y zyklus sum21
imax = 2000; %Number of iteration cycles
for ii = 1:imax    %loop over the iterations
    u(1:100) = uold(1:100) +  beta * inv(G) * eold(2:101);
    u(101) = 0;   % not needed relative degree r=1
    [y,t] = lsim(tfun_d(3),u,t,x0,'zoh');
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
figure(3)
% Plot the results of the error and the respsonse after imax iterations
subplot(425),semilogy(zyklus,sum21),hold
xlabel('$j$','interpreter','latex'),ylabel('$||e(k)||^2_2$','interpreter','latex')
strbeta = num2str(beta);
title(['Inversion Based: \beta = ' strbeta])
axis([0 200 1e-30 1e2])
subplot(426),plot(t/T,y,'g',t/T,r,'--r'),xlabel('$k$','interpreter','latex')
ylabel('- $y(k)$, - - - $r(k)$','interpreter','latex')
axis([0 100 0 5])
title('{\bf $G(s) = \frac{0.5}{(0.2 s + 1)} $}','interpreter','latex')

%% Inversion based - 4
% Initialize the plant input
uold = zeros(101,1);
u = zeros(101,1);
% Initialize the error
eold = zeros(101,1);
e = zeros(101,1);
% for a safe model beta<2
beta = 1/3.7;
% beta = 0.8;
%
clear y zyklus sum21
imax = 2000; %Number of iteration cycles
for ii = 1:imax    %loop over the iterations
    u(1:100) = uold(1:100) +  beta * inv(G) * eold(2:101);
    u(101) = 0;   % not needed relative degree r=1
    [y,t] = lsim(tfun_d(4),u,t,x0,'zoh');
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
figure(3)
% Plot the results of the error and the respsonse after imax iterations
subplot(427),semilogy(zyklus,sum21),hold
xlabel('$j$','interpreter','latex'),ylabel('$||e(k)||^2_2$','interpreter','latex')
strbeta = num2str(beta);
title(['Inversion Based: \beta = ' strbeta])
axis([0 200 1e-30 1e2])
subplot(428),plot(t/T,y,'g',t/T,r,'--r'),xlabel('$k$','interpreter','latex')
ylabel('- $y(k)$, - - - $r(k)$','interpreter','latex')
axis([0 100 0 5])
title('{\bf $G(s) = \frac{0.5}{(3 s + 1)} $}','interpreter','latex')

saveas(figure (3),'Q2_inversion_ILC','epsc')
