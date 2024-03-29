function [data, Arsig, x, lambdamax] = gen_ar(M, N, P, K, perc, Arsig)

% P=10; %order of AR-model
% N=1000; %number of data-points
% M=10; %number of channels;
sigma=.1; %scale of random AR-parameters
% K = ceil(M^2/10);

N0=1000; %length of ignored start 

if nargin < 5
  perc = 1;
end

if length(K) == 1
    inddiag = linspace(1, M^2, M);
    indndiag = setdiff(1:M^2, inddiag);
    per = randperm(length(indndiag));
    indndiag = indndiag(per);
    indndiag = indndiag(1:K);
    ind = [inddiag indndiag];
else
    inddiag = linspace(1, M^2, M);
    ind = unique([inddiag K]);
end
% RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)))

if nargin < 6
  lambdamax=10;
  while lambdamax> 1
      Arsig=[];
      for k=1:P
        aloc = zeros(M);
        aloc(ind) = double(rand(length(ind), 1) < perc).*randn(length(ind), 1)*sigma;
        Arsig=[Arsig,aloc];
      end
      E=eye(M*P);
      AA=[Arsig;E(1:end-M,:)];
      lambda=eig(AA);
      lambdamax=max(abs(lambda));
  end
end

x=randn(M,N+N0);
y=x;
for i=P+1:N+N0;
    yloc=reshape(fliplr(y(:,i-P:i-1)),[],1);
    y(:,i)=Arsig*yloc+x(:,i);
end
data=y(:,N0+1:end);

x = x(:, N0+1:end);

Arsig = reshape(Arsig, M, M, P);




