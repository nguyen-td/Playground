%% proper multiplication in matrix form
A = randn(2,4);
B = randn(4,6);

C = A * B;

%% as sums
C1 = zeros(2,6);
for i = 1:2
    for j = 1:6
        C1_sum = 0;
        for k = 1:4
             C1_sum = C1_sum + A(i,k) * B(k,j);
        end
        C1(i,j) = C1_sum;
    end
end
C(2,2) == C1(2,2);