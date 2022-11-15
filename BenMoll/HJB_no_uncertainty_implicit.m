clear all; clc;

tic;

s = 2;
r = 0.045;
rho = 0.05;
w = .1;

I=500;
amin = -0.02;
amax = 1;
a = linspace(amin,amax,I)';
da = (amax-amin)/(I-1);

maxit=10000;
crit = 10^(-6);

Delta = 1000;

dVf = zeros(I,1);
dVb = zeros(I,1);
c = zeros(I,1);

%% INITIAL GUESS
v0 = (w + r.*a).^(1-s)/(1-s)/rho;
v = v0;

for n=1:maxit
    V = v;
    % forward difference
    dVf(1:I-1) = (V(2:I)-V(1:I-1))/da;
    dVf(I) = (w + r.*amax).^(-s); %will never be used, but impose state constraint a<=amax just in case
    % backward difference
    dVb(2:I) = (V(2:I)-V(1:I-1))/da;
    dVb(1) = (w + r.*amin).^(-s); %state constraint boundary condition
    
    %consumption and savings with forward difference
    cf = dVf.^(-1/s);
    ssf = w + r.*a - cf;
    %consumption and savings with backward difference
    cb = dVb.^(-1/s);
    ssb = w + r.*a - cb;
    %consumption and derivative of value function at steady state
    c0 = w + r.*a;
    dV0 = c0.^(-s);
    
    %% dV_upwind makes a choice of forward or backward differences based on
    % the sign of the drift    
    If = ssf > 0; %positive drift --> forward difference
    Ib = ssb < 0; %negative drift --> backward difference
    I0 = (1-If-Ib); %at steady state
    
    dV_Upwind = dVf.*If + dVb.*Ib + dV0.*I0; %important to include third term
    
    c = dV_Upwind.^(-1/s);
    u = c.^(1-s)/(1-s);

    %% CONSTRUCT MATRIX
    X = -min(ssb,0)/da;
    Y = -max(ssf,0)/da + min(ssb,0)/da;
    Z = max(ssf,0)/da;
    
    %full matrix: slower
    %     for i=2:I-1
    %         A(i,i-1) = x(i);
    %         A(i,i) = y(i);
    %         A(i,i+1) = z(i);
    %     end
    %     A(1,1)=y(1); A(1,2) = z(1);
    %     A(I,I)=y(I); A(I,I-1) = x(I);
   
    %% sparse matrix: faster
    A =spdiags(Y,0,I,I)+spdiags(X(2:I),-1,I,I)+spdiags([0;Z(1:I-1)],1,I,I);
    
    if max(abs(sum(A,2)))>10^(-12)
        disp('Improper Transition Matrix')
        break
    end
    
    B = (rho + 1/Delta)*speye(I) - A;
    
    b = u + V/Delta;
    V = B\b; %SOLVE SYSTEM OF EQUATIONS
    Vchange = V - v;
    v = V;
    dist(n) = max(abs(Vchange));
    if dist(n)<crit
        disp('Value Function Converged, Iteration = ')
        disp(n)
        break
    end
end
toc;

%% Graphs
set(gca,'FontSize',14)
plot(dist,'LineWidth',2)
grid
xlabel('Iteration')
ylabel('||V^{n+1} - V^n||')

Verr = c.^(1-s)/(1-s) + dVb.*(w + r.*a - c) - rho.*V;

set(gca,'FontSize',14)
plot(a,Verr,'LineWidth',2)
grid
xlabel('k')
ylabel('Error in HJB Equation')
xlim([amin amax])

adot = w + r.*a - c;

set(gca,'FontSize',12)
plot(a,V,'LineWidth',2)
grid
xlabel('a')
ylabel('V(a)')
xlim([amin amax])

set(gca,'FontSize',14)
plot(a,c,'LineWidth',2)
grid
xlabel('a')
ylabel('c(a)')
xlim([amin amax])

set(gca,'FontSize',14)
plot(a,adot,a,zeros(1,I),'--','LineWidth',2)
grid
xlabel('a')
ylabel('s(a)')
xlim([amin amax])
