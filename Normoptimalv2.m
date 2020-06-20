%% Normoptimalv2.m
% Normoptimale ILC with a first order system 
% setup
clear e2 il u2
% System  (A=0.5; B=1; c=1; d=0;)
T=.01;
sysd = tf([1 ],[1 -.5],T);
[A,b,c,d] = ssdata(sysd);
%% Calculate the first column of G and Toeplitz
    g(1) = c*b;
    M = ones(1,1);
    for ii = 2:100
       M = M*A;
       g(ii) = c'*M*b;
    end
% Generate the  full Toeplitz Matrix
    Gvoll = toeplitz(g);
% Extract the lower trinabular part
    G = tril(Gvoll);
%% Set the weighting matrices and Learn factor L
Wdelta = eye(100);
We = eye(100);
% Learn factor
L = inv(Wdelta + G'*We*G) * G' * We;
x0 = [0];
% Reference Input and time is 101 long
r = [3*ones(1,20) 1*ones(1,20) 4 *ones(1,20) 2*ones(1,20)  1*ones(1,21) ]';
t = [0:.01:1]';
% run the iterations
% Initialize the plant input
uold = zeros(101,1);
u = zeros(101,1);
% Initialize the error
eold = zeros(101,1);
e = zeros(101,1);
% run the iteration loop
imax = 150; %Number of iteration cycles
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
subplot(221),semilogy(il,e2)
xlabel('$j$','interpreter','latex'),ylabel('$||e(k)||^2_2$','interpreter','latex')
title('${\bf W}_{\! \! \Delta u}=w_{\Delta u} {\bf I}$, $w_{\Delta u}=1$, ${\bf W}_{\! e}={\bf I}$','interpreter','latex')
axis([0 150 1e-32 1e4])
subplot(222),plot(t/T,y,'g',t/T,r,'--r'),xlabel('$k$','interpreter','latex')
ylabel('- $y(k)$, - - - $r(k)$','interpreter','latex')
axis([0 100 0 5])
%% Calculate the eigenvalues and singular values 
% as per \ref{l195} (eq 4.29) and \ref{l196} (eq 4.30) of notes
%
ME = eye(100) - G * L;
MU = eye(100) - L * G;
ME_eigmax = max( abs( eig(ME) )  )
ME_svdmax = max( abs( svd(ME) )  )
MU_eigmax = max( abs( eig(MU) )  )
MU_svdmax = max( abs( svd(MU) )  )
ME_eig(1) = ME_eigmax ;
ME_svd(1) = ME_svdmax ;
MU_svd(1) = MU_svdmax ;
%
subplot(223),semilogy(il,e2,'b'),hold
xlabel('$j$','interpreter','latex'),ylabel('$||e(k)||^2_2$','interpreter','latex')
axis([0 150 1e-32 1e4])
subplot(224),semilogx(1,ME_eig(1),'*b'),hold
xlabel('$w_{\! \Delta u}$','interpreter','latex'),ylabel('$|\lambda|_{max}$','interpreter','latex')
%% **********************************************************
% new weighting factors
%**********************************************************
Wdelta = 0.1*eye(100);
We = eye(100);
% Learn factor
L = inv(Wdelta + G'*We*G) * G' * We;
x0 = [0];
% Reference Input and time is 101 long
r = [3*ones(1,20) 1*ones(1,20) 4 *ones(1,20) 2*ones(1,20)  1*ones(1,21) ]';
t = [0:.01:1]';
% plant input
uold = zeros(101,1);
u = zeros(101,1);
% Error
eold = zeros(101,1);
e = zeros(101,1);
%
imax = 150; %Number of iteration cycles
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
%
% Calculate the eigenvalues and singular values 
% as per \ref{l195} (eq 4.29) and \ref{l196} (eq 4.30) of notes
%
ME = eye(100) - G * L;
MU = eye(100) - L * G;
ME_eigmax = max( abs( eig(ME) )  );
ME_svdmax = max( abs( svd(ME) )  );
MU_eigmax = max( abs( eig(MU) )  );
MU_svdmax = max( abs( svd(MU) )  );
ME_eig(2) = ME_eigmax ;
ME_svd(2) = ME_svdmax ;
MU_svd(2) = MU_svdmax ;
%
subplot(223),semilogy(il,e2,'g')
subplot(224),semilogx(.1,ME_eig(2),'*g')
%% **********************************************************
% new weighting factors
%**********************************************************
Wdelta = 0.01*eye(100);
We = eye(100);
% Learn factor
L = inv(Wdelta + G'*We*G) * G' * We;
x0 = [0];
% Reference Input and time is 101 long
r = [3*ones(1,20) 1*ones(1,20) 4 *ones(1,20) 2*ones(1,20)  1*ones(1,21) ]';
t = [0:.01:1]';
% plant input
uold = zeros(101,1);
u = zeros(101,1);
% Error
eold = zeros(101,1);
e = zeros(101,1);
%
imax = 150; %Number of iteration cycles
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

% Calculate the eigenvalues and singular values 
% as per \ref{l195} (eq 4.29) and \ref{l196} (eq 4.30) of notes
%
ME = eye(100) - G * L;
MU = eye(100) - L * G;
ME_eigmax = max( abs( eig(ME) )  );
ME_svdmax = max( abs( svd(ME) )  );
MU_eigmax = max( abs( eig(MU) )  );
MU_svdmax = max( abs( svd(MU) )  );
ME_eig(3) = ME_eigmax ;
ME_svd(3) = ME_svdmax ;
MU_svd(3) = MU_svdmax ;
%
subplot(223),semilogy(il,e2,'r')
subplot(224),semilogx(.01,ME_eig(3),'*r')
%% **********************************************************
% new weighting factors
%**********************************************************
Wdelta = 10*eye(100);
We = eye(100);
% Learn factor
L = inv(Wdelta + G'*We*G) * G' * We;
x0 = [0];
% Reference Input and time is 101 long
r = [3*ones(1,20) 1*ones(1,20) 4 *ones(1,20) 2*ones(1,20)  1*ones(1,21) ]';
t = [0:.01:1]';
% plant input
uold = zeros(101,1);
u = zeros(101,1);
% Error
eold = zeros(101,1);
e = zeros(101,1);
%
imax = 150; %Number of iteration cycles
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
% Calculate the eigenvalues and singular values 
% as per \ref{l195} (eq 4.29) and \ref{l196} (eq 4.30) of notes
%
ME = eye(100) - G * L;
MU = eye(100) - L * G;
ME_eigmax = max( abs( eig(ME) )  );
ME_svdmax = max( abs( svd(ME) )  );
MU_eigmax = max( abs( eig(MU) )  );
MU_svdmax = max( abs( svd(MU) )  );
ME_eig(4) = ME_eigmax ;
ME_svd(4) = ME_svdmax ;
MU_svd(4) = MU_svdmax ;
%
subplot(223),semilogy(il,e2,'k')
subplot(224),semilogx(10,ME_eig(4),'*k')
axis([0.01 10 0 1])
%print -depsc2 normoptimal.eps
%% extra
figure
semilogx([0.01 0.1 1 10],ME_eig-ME_svd,'b',[0.01 0.1 1 10],ME_eig-MU_svd,'r')
%semilogx([0.01 0.1 1 10],ME_eig,'b',[0.01 0.1 1 10],ME_svd,'g',[0.01 0.1 1 10],MU_svd,'r')

    