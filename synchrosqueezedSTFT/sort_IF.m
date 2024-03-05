function I = sort_IF(tfr,tf,Pnei,X,M,L,Ncomp,N,sumM)



tf = min(max(tf,1),M/2);
[mask] = compMask(round(tf),Pnei,M/2,sumM);
x_hat = zeros(Ncomp,N);
for c = 1:Ncomp
   x_hat(c,:) = real(rectfrgab(tfr .* mask(:,:,c), L, M));
end

% Match components and reordering for comparison
[I,~] = match_components(X, x_hat); 