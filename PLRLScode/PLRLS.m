function [Z,S,U,err] =PLRLS(X,Z_ini,lambda1,lambda2,beta,max_iter,alpha,yn)
% The code is a part of main function
%% Input X: genenum*cellnum, Ctg=inv(X'*X+2*eye(size(X,2)))
% Z_ini cellnum*cellnum
% positive paramters: lambda1,lambda2,max_iter
[m,n] = size(X);
% ---------- Initilization -------- %
miu = 0.01;
rho = 1.2;
max_miu = 1e8;
tol  = 1e-5;
tol2 = 1e-2;
C1 = zeros(n,n);
C2 = zeros(n,n);
%  H= zeros(n,n);%ILRR（Without H,beta=0）
H= fraction(X,0,yn);
% H= fraction_L(abs(X),0,yn);
% H= cosSim(X');
% H=L_chebyshev(X) ;
%第二个参数(0--functional functin;2--euclidean;others--chebychev)
%第三个参数决定分数函数的k=0.25/0.5/0.75
% H=corr(X,'type','Pearson');
% H=simlarityGS(X);%

distX= fraction(X,0,0.25);
% distX = L2_distance_1(X,X);%欧氏距离的平方
distX=(1-alpha)*distX/max(max(distX))+alpha*cosSim(X');
% distX=zeros(n,n);
for iter = 1:max_iter
    if iter == 1
        Z = Z_ini;
        S = Z_ini;
        U = Z_ini;
        clear Z_ini F_ini
    end
    S_old = S;
    U_old = U;
    Z_old = Z;
    % -------- Update Z --------- %
    err(iter)=norm(X-X*Z,'fro');
    T=lambda2*(X'*X);
    Z=(beta*eye(n)+T+2*miu*eye(n))\(beta*H+T+miu*(U+S)-C1-C2);
    Z = Z- diag(diag(Z));
    % -------- Update S --------- %
    dist  = distX;
    S     = Z+(C2-dist)/miu;
    S     = S - diag(diag(S));
    for ic = 1:n
        idx    = 1:n;
        idx(ic) = [];
        S(ic,idx) = EProjSimplex_new(S(ic,idx));          %
    end
    
    % -------- Update U --------- %
    [AU,SU,VU] = svd(Z+C1/miu,'econ');
    AU(isnan(AU)) = 0;
    VU(isnan(VU)) = 0;
    SU(isnan(SU)) = 0;
    SU = diag(SU);
    SVP = length(find(SU>lambda1/miu));
    if SVP >= 1
        SU = SU(1:SVP)-lambda1/miu;
    else
        SVP = 1;
        SU = 0;
    end
    U = AU(:,1:SVP)*diag(SU)*VU(:,1:SVP)';
    %      U=(miu*Z+C3)/(2*lambda1+miu);  %% F norm
    % -------- Update C1 C2 miu -------- %
    L1 = Z-U;
    L2 = Z-S;
    C1 = C1+miu*L1;%----U
    C2 = C2+miu*L2;
    
    LL1 = norm(Z-Z_old,'fro');%返回矩阵(Z-Z_old)的 Frobenius 范数
    LL2 = norm(S-S_old,'fro');
    LL3 = norm(U-U_old,'fro');
    SLSL = max(max(LL1,LL2),LL3)/norm(X,'fro');
    if miu*SLSL < tol2
        miu = min(rho*miu,max_miu);
    end
    % --------- obj ---------- %
    leq1 = max(max(abs(L1(:))),max(abs(L2(:))));
%     stopC = max(leq1,max(abs(L3(:))));
    
    if leq1 < tol
        iter
        break;
    end
    
end
end


function d = L2_distance_1(a,b)
% compute squared Euclidean distance欧氏距离的平方
% ||A-B||^2 = ||A||^2 + ||B||^2 - 2*A'*B
% a,b: two matrices. each column is a data
% d:   distance matrix of a and b
if (size(a,1) == 1)
    a = [a; zeros(1,size(a,2))];
    b = [b; zeros(1,size(b,2))];
end
aa=sum(a.*a); bb=sum(b.*b); ab=a'*b;
d = repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab;
d = real(d);
d = max(d,0);
% % force 0 on the diagonal?
% if (df==1)
%   d = d.*(1-eye(size(d)));
% end
end

function result=cosSim(data)
%COSSIM Summary of this function goes here余弦相似度
%Detailed explanation goes here
rows=size(data,1);
result=zeros(rows,rows);
for i=1:rows
    for j=1:i
        if (norm(data(i,:))*norm(data(j,:))==0)
            result(i,j)=0;
        else
            result(i,j)=dot(data(i,:),data(j,:))/(norm(data(i,:))*norm(data(j,:))); 
        end
        result(j,i)=result(i,j);
    end
end
end

function [x ft] = EProjSimplex_new(v, k)
%% Optimization Problem
%
%  min  1/2 || x - v||^2
%  s.t. x>=0, 1'x=k
%
if nargin < 2
    k = 1;
end
ft=1;
n = length(v);
v0 = v-mean(v) + k/n;
%vmax = max(v0);
vmin = min(v0);
if vmin < 0
    f = 1;
    lambda_m = 0;
    while abs(f) > 10^-10
        v1 = v0 - lambda_m;
        posidx = v1>0;
        npos = sum(posidx);
        g = -npos;
        f = sum(v1(posidx)) - k;
        lambda_m = lambda_m - f/g;
        ft=ft+1;
        if ft > 100
            x = max(v1,0);
            break;
        end
    end
    x = max(v1,0);
else
    x = v0;
end
end


function [E] = solve_l1l2(W,lambda)
n = size(W,2);
E = W;
for i=1:n
    E(:,i) = solve_l2(W(:,i),lambda);
end
end

function [x] = solve_l2(w,lambda)
% min lambda |x|_2 + |x-w|_2^2
nw = norm(w);
if nw>lambda
    x = (nw-lambda)*w/nw;
else
    x = zeros(length(w),1);
end
end



function W = Network_Diffusion(A, K)
%K = min(2*K, round(length(A)/10));
A = A-diag(diag(A));
P = (dominateset(double(abs(A)),min(K,length(A)-1))).*sign(A);
DD = sum(abs(P'));
P = P + (eye(length(P))+diag(sum(abs(P'))));
P = (TransitionFields(P));
[U,D] = eig(P);
d = real((diag(D))+eps);
alpha = 0.8;
beta = 2;
d = (1-alpha)*d./(1-alpha*d.^beta);

D = diag(real(d));
W = U*D*U';

W = (W.*(1-eye(length(W))))./repmat(1-diag(W),1,length(W));
D=sparse(1:length(DD),1:length(DD),DD);
W=D*(W);
W = (W+W')/2;

end



