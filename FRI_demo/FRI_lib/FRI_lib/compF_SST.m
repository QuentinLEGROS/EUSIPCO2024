function [F_sst]=compF_SST(M,mean)
%
% Compute the impulse response function (IRF) of a synchrosqueezed signal. The Gabor
% kernel is approximated using a Gaussian distribution.
%
% INPUT:
% M             : Number of frequential bin
% mean          : mean of the IRF. can be modified arbitrarily as long as
%                 it remains far from the signal edges
% 
% 
% OUTPUT:
% F_mat         : Postulated obervation model


F_sst = normpdf(-(M/4)+1:((3*M)/4),mean,0.5);
F_sst = F_sst((M/4)+1:((3*M)/4));
F_sst = (F_sst./sum(F_sst))';% Normalization