%% owens07a_no_I_term.m 
% Example from Owens07a with a differnt Model family no I-Term, to
% motivate a gradient method
% set up nominal plant
s = tf('s');
G0 = 1/(s^2);
T = 0.1;
G0d = c2d(G0,T);
%% Alternative Method - use uncertain parameters to generate a family of models
K = ureal('K',1,'range',[0.5 1.5]);
T1 = ureal('T1',1.6,'range',[0.2 3.]);
tfunsicher = tf(K,[T1 1]);
%% Now only use the four extreme case models and plot Nyquist
% for plant inversion method
% as in Harte
tfun(1) = tf(1.5,[0.2 1]);
tfun(2) = tf(1.5,[3.0 1]);
tfun(3) = tf(0.5,[0.2 1]);
tfun(4) = tf(0.5,[3.0 1]);
% convert to discrete time using ZOH
tfun_d(1) = c2d( tfun(1) , T);
tfun_d(2) = c2d( tfun(2) , T);
tfun_d(3) = c2d( tfun(3) , T);
tfun_d(4) = c2d( tfun(4) , T);
%
subplot(221),nyquist(tfun_d(1)/G0d,'b',tfun_d(2)/G0d,'r',tfun_d(3)/G0d,'g',tfun_d(4)/G0d,'k'),hold
% contained in a circle at  1/\beta with radius \beta
betainv = 3.7; %Radius
beta = 1/betainv;
betaharte = beta;
einre = [0:.001:betainv]+betainv;
lengthb = length(einre);
einim = sqrt(betainv^2 - (einre-betainv).^2);
plot(einre,einim,'--','LineWidth',2),plot(einre,-einim,'--','LineWidth',2)
plot(einre-betainv,einim(lengthb:-1:1),'--','LineWidth',2)
plot(einre-betainv,-einim(lengthb:-1:1),'--','LineWidth',2)
%axis('square')
axis([-0.3 11 -4 4])
xlabel('$Re(U)$','interpreter','latex'),ylabel('$Im(U)$','interpreter','latex')
title('Plant Inversion Method')
%% Now use the gradient method for the four extreme case models
% Calculation follows Owens
subplot(222),nyquist(G0d'*tfun_d(1),'b',G0d'*tfun_d(2),'r',G0d'*tfun_d(3),'g',G0d'*tfun_d(4),'k'),hold
betainv = 0.75; %Radius
beta = 1/betainv;
betaowens = beta;
einre = [0:.001:betainv]+betainv;
lengthb = length(einre);
einim = sqrt(betainv^2 - (einre-betainv).^2);
plot(einre,einim,'--','LineWidth',2),plot(einre,-einim,'--','LineWidth',2)
plot(einre-betainv,einim(lengthb:-1:1),'--','LineWidth',2)
plot(einre-betainv,-einim(lengthb:-1:1),'--','LineWidth',2)
xlabel('$Re(U)$','interpreter','latex'),ylabel('$Im(U)$','interpreter','latex')
%axis('square')
axis([-0.3 2 -.8 .8])
title('Gradient Method')
%% ILC Setup
% Time dependence of the 2-norm of the 4 extreme case models
% first the super vector format for the nominal model G0d is generated
% note all 4 discrete time real models have  a relative degree of 1
[a,b,c,d] = ssdata(G0d);
% Markov-Parameter of the nominal model
% The first row for the 100 time points
    g(1) = c*b;
    M = eye(1,1);
    for i = 2:100
       M = M*a;
       g(i) = c*M*b;
    end
% The full Toeplitz matrix is
    Gvoll = toeplitz(g);
% the extract the lower triangular part
    G = tril(Gvoll);
% all vectors are 101 long for a relative degree of r=1
% this cause the shift
% The reference is 101 long, the first value is unimportant
r = [1*ones(1,101) ]';
r = [3*ones(1,20) 1*ones(1,20) 4 *ones(1,20) 2*ones(1,20)  1*ones(1,21) ]';
t = [0:.1:10]';
%% ILC for the Gradient method 
% ILC - Harte-method
beta = betaharte;
x0 = [0; 0];
%% 1. real Model with ILC
% plant input
uold = zeros(101,1);
u = zeros(101,1);
% Error
eold = zeros(101,1);
%ealt = rand(101,1);
e = zeros(101,1);
% Learn factor beta from above
%
imax = 500; %Number of cycles
for i = 1:imax    %over all cycles
    u(1:100) = uold(1:100) +  beta * inv(G) * eold(2:101);
    u(101) = 0;  %not used as the realtive degree is r=1
    [y,t] = lsim(tfun_d(1),u,t,x0,'zoh');
    sum = 0;
    for k = 2:101
        e(k) = r(k) - y(k);
        sum = sum + e(k)^2;
    end
    sum21(i) = sum;
    zyklus(i) = i;
    eold = e;
    uold = u;
end
%% 2. real Model with ILC
%  plant input
uold = zeros(101,1);
u = zeros(101,1);
% Error
eold = zeros(101,1);
%ealt = rand(101,1);
e = zeros(101,1);
% Learn factor beta from above
%
imax = 500; %Number of cycles
for i = 1:imax    %over all cycles
    u(1:100) = uold(1:100) +  beta * inv(G) * eold(2:101);
    u(101) = 0;  %not used as the realtive degree is r=1
    [y,t] = lsim(tfun_d(2),u,t,x0,'zoh');
    sum = 0;
    for k = 2:101
        e(k) = r(k) - y(k);
        sum = sum + e(k)^2;
    end
    sum22(i) = sum;
    zyklus(i) = i;
    eold = e;
    uold = u;
end
%% 3. real Model with ILC
%  plant input
uold = zeros(101,1);
u = zeros(101,1);
% Error
eold = zeros(101,1);
%ealt = rand(101,1);
e = zeros(101,1);
% Learn factor beta from above
%
imax = 500; %Number of cycles
for i = 1:imax    %over all cycles
    u(1:100) = uold(1:100) + beta * inv(G) * eold(2:101);
    u(101) = 0;  %not used as the realtive degree is r=1
    [y,t] = lsim(tfun_d(3),u,t,x0,'zoh');
    sum = 0;
    for k = 2:101
        e(k) = r(k) - y(k);
        sum = sum + e(k)^2;
    end
    sum23(i) = sum;
    zyklus(i) = i;
    eold = e;
    uold = u;
end
%
% this result is shown below
u3 = u;
y3 = y;
%
%
%% 4. real Model with ILC
%  plant input
uold = zeros(101,1);
u = zeros(101,1);
% Error
eold = zeros(101,1);
%ealt = rand(101,1);
e = zeros(101,1);
% Learn factor beta from above
%
imax = 500; %Number of cycles
for i = 1:imax    %over all cycles
    u(1:100) = uold(1:100) + beta * inv(G) * eold(2:101);
    u(101) = 0;  %not used as the realtive degree is r=1
    [y,t] = lsim(tfun_d(4),u,t,x0,'zoh');
    sum = 0;
    for k = 2:101
        e(k) = r(k) - y(k);
        sum = sum + e(k)^2;
    end
    sum24(i) = sum;
    zyklus(i) = i;
    eold = e;
    uold = u;
end
%
subplot(223),semilogy(zyklus,sum21,'b',zyklus,sum22,'r',zyklus,sum23,'g',zyklus,sum24,'k')
xlabel('$j$','interpreter','latex')
ylabel('$||e(k)|| _2^2$','interpreter','latex')
axis([0 500 1e-30 10000])
title('Plant Inversion Method')
%
%
%% ILC Gradient Method setup
% use Ownes-method
%
beta = betaowens;
x0 = [0; 0];
% all vectors are 101 long for a relative degree of r=1
% this cause the shift
% The reference is 101 long, the first value is unimportant
r = [1*ones(1,101) ]';
r = [3*ones(1,20) 1*ones(1,20) 4 *ones(1,20) 2*ones(1,20)  1*ones(1,21) ]';
t = [0:.1:10]';
%
%% 1. real Model with ILC
%  plant input
uold = zeros(101,1);
u = zeros(101,1);
% Error
eold = zeros(101,1);
%ealt = rand(101,1);
e = zeros(101,1);
% Learn factor beta from above
%
imax = 500; %Number of cycles
for i = 1:imax    %over all cycles
    u(1:100) = uold(1:100) +  beta * G' * eold(2:101);
    u(101) = 0;  %not used as the realtive degree is r=1
    [y,t] = lsim(tfun_d(1),u,t,x0,'zoh');
    sum = 0;
    for k = 2:101
        e(k) = r(k) - y(k);
        sum = sum + e(k)^2;
    end
    sum21(i) = sum;
    zyklus(i) = i;
    eold = e;
    uold = u;
end
y1 = y;
%% 2. real Model with ILC
%  plant input
uold = zeros(101,1);
u = zeros(101,1);
% Error
eold = zeros(101,1);
%ealt = rand(101,1);
e = zeros(101,1);
% Learn factor beta from above
%
imax = 500; %Number of cycles
for i = 1:imax    %over all cycles
    u(1:100) = uold(1:100) +  beta * G' * eold(2:101);
    u(101) = 0;  %not used as the realtive degree is r=1
    [y,t] = lsim(tfun_d(2),u,t,x0,'zoh');
    sum = 0;
    for k = 2:101
        e(k) = r(k) - y(k);
        sum = sum + e(k)^2;
    end
    sum22(i) = sum;
    zyklus(i) = i;
    eold = e;
    uold = u;
end
y2 = y;
%% 3. real Model with ILC
%  plant input
uold = zeros(101,1);
u = zeros(101,1);
% Error
eold = zeros(101,1);
%ealt = rand(101,1);
e = zeros(101,1);
% Learn factor beta from above
%
imax = 500; %Number of cycles
for i = 1:imax    %over all cycles
    u(1:100) = uold(1:100) + beta * G' * eold(2:101);
    u(101) = 0;  %not used as the realtive degree is r=1
    [y,t] = lsim(tfun_d(3),u,t,x0,'zoh');
    sum = 0;
    for k = 2:101
        e(k) = r(k) - y(k);
        sum = sum + e(k)^2;
    end
    sum23(i) = sum;
    zyklus(i) = i;
    eold = e;
    uold = u;
end
% the results are shown below
u3 = u;
y3 = y;
%
%% 4. real Model with ILC
%  plant input
uold = zeros(101,1);
u = zeros(101,1);
% Error
eold = zeros(101,1);
%ealt = rand(101,1);
e = zeros(101,1);
% Learn factor beta from above
%
imax = 500; %Number of cycles
for i = 1:imax    %over all cycles
    u(1:100) = uold(1:100) + beta * G' * eold(2:101);
    u(101) = 0;  %not used as the realtive degree is r=1
    [y,t] = lsim(tfun_d(4),u,t,x0,'zoh');
    sum = 0;
    for k = 2:101
        e(k) = r(k) - y(k);
        sum = sum + e(k)^2;
    end
    sum24(i) = sum;
    zyklus(i) = i;
    eold = e;
    uold = u;
end
y4 = y;
%% Plot the gradient method results
subplot(224),semilogy(zyklus,sum21,'b',zyklus,sum22,'r',zyklus,sum23,'g',zyklus,sum24,'k')
xlabel('$j$','interpreter','latex')
ylabel('$||e(k)|| _2^2$','interpreter','latex')
axis([0 500 1e-30 10000])
title('Gradient Method')
% print -depsc2 Owens07a_no_I_term.eps
%
%% Plot the worst time response
%
figure
%subplot(221)
set(0,'DefaultLineLineWidth',1.0)
bk_font=16; % for a 16 point font
axes('Fontsize',bk_font, 'FontName', 'Helvetica')
plot(t/T,y4,'k',t/T,r,'--r'),xlabel('$k$','interpreter','latex')
ylabel('- y(k), - - - r(k)','interpreter','latex')
axis([0 100 0 5])
title('$G_p=0.5/(3.0 s+1)$','interpreter','latex')
%print -depsc2 Owens07a_no_I_part_time.eps