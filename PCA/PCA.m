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
%  ��beta�滻alf

[~,n] = size(X);
% Laplacian
d = sum(W);
% E=X-VQ';
% M=X-E-C/mu;%%%%2,1�������
% mu=3;
K = X'*X;
D = diag(d);
L = D - W;
I = eye(n,n);% ��λ����  ֻ�жԽ�Ԫ��Ϊ1
% [m,n]=size(X);
% C = zeros( m, n);

sk = max(eig(K));
sl = max(eig(L));%  ���������ֵ
% KL = (1 - beta)*(I - K/sk) +  beta*L/sl;
e = ones(n,1);
% KL = (1 - beta)*(I - K/sk) +  beta*(L/sl + (e*e')/n); %��14��
KL = (1 - beta)*(I - K/sk) ; %��14��
for ii=1:50
%============================ Compute the optimal V
[U0,V0] = eig(KL);%%%%%%%%  V0����ֵ���ɵĶԽ���U0�����������ɵ�������
val = diag(V0);%%%%%%%%%%������ֵ�ԽǾ���V0�ų�������
[~,ind] = sort(val);%%%%%%%��val��������Ԫ�ص�������˳��
V = U0(:,ind(1:k));%%%%������δ����k��embedding vector V (size:n*k,k:reduction dim);

% ================================Compute the optimal U
U = X*V;
% Compute the optimal C
% C=C+mu(E-X+UQ');
% % Compute the optimal mu
% mu=rho*mu;
err(ii) = norm(X-U*V','fro');
end
