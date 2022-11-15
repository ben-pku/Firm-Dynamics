%% Solve the steady state of BEMV firm dynamics
% Zhuokai Huang   edited on 2022/0408
clear
close all
clc
disp('solving steady state...')
tic 
%% input data
run ExogenousParameters.m
load Input/EstimatedParams.mat
Params.zbar = 0;
    % Options
options.Hopenhayn           = 0;
options.HopenhaynNumerical  = 0;
options.Transition          = 0;
%% solve steady state
[~,~,~, SUntrunc_znmat, SnUntrunc_znmat, S_znmat, Sn_znmat, v_znmat, ...
    g_znmat, Gn_znmat, Gv_znmat, q, phi, p, u, L, pi0Trunc_zn, NumGrids, Derivatives]...
    = SolveBEMV( Params, ExogParams, options ) ;

%% plot
% S_znmat 
f1 =figure(1);
surf( NumGrids.z, NumGrids.n, S_znmat);
xlabel('Productivity log(z)','FontSize',14);
ylabel('Size log(n)' ,'FontSize',14);
zlabel('Joint Surplus S','FontSize',14 );
print(f1, './surplus_steadystate','-depsc');
% g_znmat
f2 = figure(2);
surf(NumGrids.z, NumGrids.n, g_znmat);
xlabel('Productivity log(z)','FontSize',14);
ylabel('Size log(n)','FontSize',14);
zlabel('Density g(z,n)','FontSize',14 );
print(f2, './density_steadystate','-depsc');