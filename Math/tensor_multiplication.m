%% tprod documentation
% X = randn(10,9,8,7,6);
% Z = tprod(X,[1 2 -3 -4 3],X,[1 2 -3 -4 3]); % accumulate away dimensions 3&4 and squeeze the result to 3d
% MATLAB equivalent; for i=1:size(X,5); Xi=reshape(X(:,:,:,:,i),[10*9 8*7]); Z(:,:,i) = reshape(sum(Xi.*Xi,2),[10 9]); end;

%% tensorproduct
X = randn(10,9,8,7,6);
Z = zeros(10,9,6);
for i=1:size(X,5) % 6
    Xi = reshape(X(:,:,:,:,i),[10*9 8*7]);
    Z(:,:,i) = reshape(sum(Xi.*Xi,2),[10 9]);
end

%% equivalent to the thing above
H = zeros(10,9,6);
for i = 1:10
    for j = 1:9
        for m = 1:6
            for n = 1:8
                for p = 1:7
                    H(i,j,m) = X(i,j,n,p,m) * X(i,j,n,p,m);
                end
            end
        end
    end
end

h_11 = sum(H(1,1,:))
z_11 = sum(Z(1,1,:))

%% Sanity check
p = 10;
q = 12;
r = 15;
bs_est_pqr = 0;
for i = 1:n
    for j = 1:n
        for k = 1:n
            bs_est_pqr = bs_est_pqr + (a(p,i) * a(q,j) * a(r,k) * d(i,j,k));
        end
    end
end
disp(bs_est_pqr)

% tensor notation
bs_est = tensorprod(a, tensorprod(a, tensorprod(a, d, 2, 3), 2, 3), 2, 3);
figure; imagesc(squeeze(abs(bs_orig(:,:,5))))
figure; imagesc(squeeze(abs(bs_est(:,:,5))))
disp(bs_est(p,q,r))