function [F_mat] = compute_F(M,L)
% [F_mat] = compute_F(M,L)
%
% Compute the data distribution and its shifted version
%
% 
% INPUT:
% M      : number of frequency bins to process
% L      : window time spread parameter
%
% OUTPUT:
% F_mat    : convolution matrix of the data distribution
%
% Author: Q.Legros (quentin.legros@telecom-paris.fr)
% Date: 12-may-2021


val = transpose(Fh(-(M/4)+1:((3*M)/4), M, L ));

F_mat = zeros(M,M/2);
F_mat(:,1) = val;
for i = 1:(M/2)-1
    F_mat(:,i+1) = [F_mat(end,i);F_mat(1:end-1,i)];
end
% F_mat = F_mat((M/4)+1:((3*M)/4),:); % Truncation to the same lenght than the data
F_mat = F_mat((M/4):((3*M)/4)-1,:); % Truncation to the same lenght than the data
F_mat = F_mat./max(F_mat);% Normalization
