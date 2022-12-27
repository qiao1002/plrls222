function [U,V,err] = GLPCA(X,W,beta,k)
%% Graph Laplacian PCA function by Bo Jiang in Anhui university  2015/7/8
% input: vector data X (size:p*n,p:orignal data dim,n:data num); 
%        graph data W(size:n*n)
%        beta:weighting parameter
%        k:reduced dimension
% output:embedding vector V (size:n*k,k:reduction dim);
%        projection U (size:p*k)
% obj: J = \|X-UV^T\|^2_F+alf*trV^T*(D-W)*V    s.t. V^T*V=I
% note that the parameter beta is the replacement of balanced parameter alf, see Eq.(10) in our paper
%  用beta替换alf

[~,n] = size(X);
% Laplacian
d = sum(W);
% E=X-VQ';
% M=X-E-C/mu;%%%%2,1范数求解
% mu=3;
K = X'*X;
D = diag(d);
L = D - W;
I = eye(n,n);% 单位矩阵  只有对角元素为1
% [m,n]=size(X);
% C = zeros( m, n);

sk = max(eig(K));
sl = max(eig(L));%  求最大特征值
% KL = (1 - beta)*(I - K/sk) +  beta*L/sl;
e = ones(n,1);
% KL = (1 - beta)*(I - K/sk) +  beta*(L/sl + (e*e')/n); %（14）
KL = (1 - beta)*(I - K/sk) ; %（14）
for ii=1:50
%============================ Compute the optimal V
[U0,V0] = eig(KL);%%%%%%%%  V0特征值构成的对角阵，U0特征向量构成的列向量
val = diag(V0);%%%%%%%%%%把特征值对角矩阵V0排成列向量
[~,ind] = sort(val);%%%%%%%求val列向量里元素的索引，顺序
V = U0(:,ind(1:k));%%%%这里面未定义k；embedding vector V (size:n*k,k:reduction dim);

% ================================Compute the optimal U
U = X*V;
% Compute the optimal C
% C=C+mu(E-X+UQ');
% % Compute the optimal mu
% mu=rho*mu;
err(ii) = norm(X-U*V','fro');
end
