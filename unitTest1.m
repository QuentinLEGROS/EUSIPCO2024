clear 
close all

addpath('/home/dfourer/Recherche/Publications/30_EUSIPCO_2017/ASTRES_toolbox/RLMR/');


N = 500;
M = 500;
Ncomp = 1;
M0 = 30;

k = 3;
n_pad = 50;
L = 20;
[Fr,a,b] = init_recursif_data(M,L,k);

x = randn(1, N).';
[tf, tfrp] = RecursiveFRI(x,k,L,n_pad,M,N,a,b,Ncomp,Fr,M0);

[tfr0, nfreqs] = recursive_stft(x, k, L, 0, M/2-1, M);