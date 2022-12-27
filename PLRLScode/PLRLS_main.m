clear
clc
load('Data_Treutlein.mat')%1.1,0,-0.2
% load('Data_Ting.mat')%0,0,1
lambda=10^(1);
gamma=10^(2);
[in_X,] = FilterGenesZero(in_X);
fea=in_X;
gnd=true_labs(:);
num_class =length(unique(gnd)); % The number of classes
% ---------- Initilization for Z  -------- %
fea = double(fea);  %% 
select_sample = [];
select_gnd    = [];
for i = 1:num_class  %% 
    idx = find(gnd == i);%把gnd=i的索引放在idx中
    idx_sample    = fea(idx,:);%在原始矩阵中找到类别为i的所有细胞
    select_sample = [select_sample;idx_sample];
    select_gnd    = [select_gnd;gnd(idx)];
end
fea = select_sample';%select_sample，select_gnd已经按照类别从小到大排好序了。
fea = fea./repmat(sqrt(sum(fea.^2)),[size(fea,1) 1]);%归一化
gnd = select_gnd;
X = fea;  %% 
clear fea select_gnd select_sample idx
%%  ---------- Initilization for Z  -------- %
options = [];
options.NeighborMode = 'KNN';
options.k = 10;
%options.WeightMode = 'Binary';      % Binary  HeatKernel
options.WeightMode = 'Cosine';
Z = constructW(X',options);
Z = full(Z);
Z1 = Z-diag(diag(Z));
Z_ini = (Z1+Z1')/2;
clear LZ DZ Z fea Z1
H= fraction(X,0,0.25);
%%  ---------- Optimization for Z  -------- %
max_iter= 80; %default value
% lambda=10^(1.4);
% gamma=10^(1);
% beta=10^(0);%-5-0
% alpha=0.7;
ii=0;
for jj=(1.1):1:(1.1)
    beta=10^(0);
    ii=ii+1;
    i=0;
    for j=(-0.2):0.1:(-0.2)
       i=i+1;
       alpha=0.7;
       [Z,S,U,err] = ISLRR(X,Z_ini,lambda,gamma,beta,max_iter,alpha,0.25);
        similarity=(abs(Z)+abs(Z'))/2;
%         NUMC = 2 : 15;
%         [K1, ~, K12, ~] = Estimate_Number_of_Clusters_given_graph(similarity, NUMC);
%         fprintf('The number of clusters estimated by eigen gap : %d\n', K1);
%         fprintf('The number of clusters estimated by rotation cost : %d\n', K12);
        [result_label, kerNS]= SpectralClustering(similarity,num_class);
        % ---------- Evaluation measurement  -------- %
        NMI=Cal_NMI_newused(gnd, result_label);
        ARI=Contingency_ARI_newused(gnd, result_label);
        [result1(ii,i)]=NMI;
        [result2(ii,i)]=ARI;
    end
end
result1
[maxr,index]=max(result1)
[NMI,Purity,ARI]
%% 
function [groups, kerNS] = SpectralClustering(CKSym,n)
%--------------------------------------------------------------------------
% This function takes an adjacency matrix of a graph and computes the
% clustering of the nodes using the spectral clustering algorithm of
% Ng, Jordan and Weiss.
% CMat: NxN adjacency matrix
% n: number of groups for clustering
% groups: N-dimensional vector containing the memberships of the N points
% to the n groups obtained by spectral clustering
%--------------------------------------------------------------------------
% Copyright @ Ehsan Elhamifar, 2012
% Modified @ Chong You, 2015
%--------------------------------------------------------------------------
warning off;
N = size(CKSym,1);
% Normalized spectral clustering according to Ng & Jordan & Weiss
% using Normalized Symmetric Laplacian L = I - D^{-1/2} W D^{-1/2}
DN = diag(1./sqrt(sum(CKSym)+eps) );   %eps=2.2204e-16
LapN = speye(N) - DN * CKSym * DN;  % speye(N)生成N*N对角元素为1的矩阵  构建拉普拉斯矩阵
[~,~,vN] = svd(LapN); %% 奇异值分解
kerN = vN(:,N-n+1:N);
%kerN = vN(:,N-12:N);
normN = sum(kerN .^2, 2) .^.5;%% normalize the matrix U by L2-Norm
kerNS = bsxfun(@rdivide, kerN, normN + eps);  %% 矩阵kerN 每一行除以normN+eps
%-------------
%Y = pdist(kerNS);
%Z = linkage(Y);
%groups = cluster(Z,'maxclust',n);
MAXiter = 1000; % Maximum number of iterations for KMeans
REPlic = 100; % Number of replications for KMeans
groups = kmeans(kerNS,n,'maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');  %% 调用kmeans聚类
end

function W = constructW(fea,options)
%	Usage:
%	W = constructW(fea,options)
%
%	fea: Rows of vectors of data points. Each row is x_i
%   options: Struct value in Matlab. The fields in options that can be set:
%
%           NeighborMode -  Indicates how to construct the graph. Choices
%                           are: [Default 'KNN']
%                'KNN'            -  k = 0
%                                       Complete graph
%                                    k > 0
%                                      Put an edge between two nodes if and
%                                      only if they are among the k nearst
%                                      neighbors of each other. You are
%                                      required to provide the parameter k in
%                                      the options. Default k=5.
%               'Supervised'      -  k = 0
%                                       Put an edge between two nodes if and
%                                       only if they belong to same class.
%                                    k > 0
%                                       Put an edge between two nodes if
%                                       they belong to same class and they
%                                       are among the k nearst neighbors of
%                                       each other.
%                                    Default: k=0
%                                   You are required to provide the label
%                                   information gnd in the options.
%
%           WeightMode   -  Indicates how to assign weights for each edge
%                           in the graph. Choices are:
%               'Binary'       - 0-1 weighting. Every edge receiveds weight
%                                of 1.
%               'HeatKernel'   - If nodes i and j are connected, put weight
%                                W_ij = exp(-norm(x_i - x_j)/2t^2). You are
%                                required to provide the parameter t. [Default One]
%               'Cosine'       - If nodes i and j are connected, put weight
%                                cosine(x_i,x_j).
%
%            k         -   The parameter needed under 'KNN' NeighborMode.
%                          Default will be 5.
%            gnd       -   The parameter needed under 'Supervised'
%                          NeighborMode.  Colunm vector of the label
%                          information for each data point.
%            bLDA      -   0 or 1. Only effective under 'Supervised'
%                          NeighborMode. If 1, the graph will be constructed
%                          to make LPP exactly same as LDA. Default will be
%                          0.
%            t         -   The parameter needed under 'HeatKernel'
%                          WeightMode. Default will be 1
%         bNormalized  -   0 or 1. Only effective under 'Cosine' WeightMode.
%                          Indicates whether the fea are already be
%                          normalized to 1. Default will be 0
%      bSelfConnected  -   0 or 1. Indicates whether W(i,i) == 1. Default 0
%                          if 'Supervised' NeighborMode & bLDA == 1,
%                          bSelfConnected will always be 1. Default 0.
%            bTrueKNN  -   0 or 1. If 1, will construct a truly kNN graph
%                          (Not symmetric!). Default will be 0. Only valid
%                          for 'KNN' NeighborMode
%
%
%    Examples:
%
%       fea = rand(50,15);
%       options = [];
%       options.NeighborMode = 'KNN';
%       options.k = 5;
%       options.WeightMode = 'HeatKernel';
%       options.t = 1;
%       W = constructW(fea,options);
%
%
%       fea = rand(50,15);
%       gnd = [ones(10,1);ones(15,1)*2;ones(10,1)*3;ones(15,1)*4];
%       options = [];
%       options.NeighborMode = 'Supervised';
%       options.gnd = gnd;
%       options.WeightMode = 'HeatKernel';
%       options.t = 1;
%       W = constructW(fea,options);
%
%
%       fea = rand(50,15);
%       gnd = [ones(10,1);ones(15,1)*2;ones(10,1)*3;ones(15,1)*4];
%       options = [];
%       options.NeighborMode = 'Supervised';
%       options.gnd = gnd;
%       options.bLDA = 1;
%       W = constructW(fea,options);
%
%
%    For more details about the different ways to construct the W, please
%    refer:
%       Deng Cai, Xiaofei He and Jiawei Han, "Document Clustering Using
%       Locality Preserving Indexing" IEEE TKDE, Dec. 2005.
%
%
%    Written by Deng Cai (dengcai2 AT cs.uiuc.edu), April/2004, Feb/2006,
%                                             May/2007
%
bSpeed  = 1;
if (~exist('options','var'))
    options = [];
end
if isfield(options,'Metric')
    warning('This function has been changed and the Metric is no longer be supported');
end
if ~isfield(options,'bNormalized')
    options.bNormalized = 0;
end
%=================================================
if ~isfield(options,'NeighborMode')
    options.NeighborMode = 'KNN';
end
switch lower(options.NeighborMode)
    case {lower('KNN')}  %For simplicity, we include the data point itself in the kNN
        if ~isfield(options,'k')
            options.k = 5;
        end
    case {lower('Supervised')}
        if ~isfield(options,'bLDA')
            options.bLDA = 0;
        end
        if options.bLDA
            options.bSelfConnected = 1;
        end
        if ~isfield(options,'k')
            options.k = 0;
        end
        if ~isfield(options,'gnd')
            error('Label(gnd) should be provided under ''Supervised'' NeighborMode!');
        end
        if ~isempty(fea) && length(options.gnd) ~= size(fea,1)
            error('gnd doesn''t match with fea!');
        end
    otherwise
        error('NeighborMode does not exist!');
end
%=================================================
if ~isfield(options,'WeightMode')
    options.WeightMode = 'HeatKernel';
end
bBinary = 0;
bCosine = 0;
switch lower(options.WeightMode)
    case {lower('Binary')}
        bBinary = 1;
    case {lower('HeatKernel')}
        if ~isfield(options,'t')
            nSmp = size(fea,1);
            if nSmp > 3000
                D = EuDist2(fea(randsample(nSmp,3000),:));
            else
                D = EuDist2(fea);
            end
            options.t = mean(mean(D));
        end
    case {lower('Cosine')}
        bCosine = 1;
    otherwise
        error('WeightMode does not exist!');
end
%=================================================
if ~isfield(options,'bSelfConnected')
    options.bSelfConnected = 0;
end
%=================================================
if isfield(options,'gnd')
    nSmp = length(options.gnd);
else
    nSmp = size(fea,1);
end
maxM = 62500000; %500M
BlockSize = floor(maxM/(nSmp*3));


if strcmpi(options.NeighborMode,'Supervised')
    Label = unique(options.gnd);
    nLabel = length(Label);
    if options.bLDA
        G = zeros(nSmp,nSmp);
        for idx=1:nLabel
            classIdx = options.gnd==Label(idx);
            G(classIdx,classIdx) = 1/sum(classIdx);
        end
        W = sparse(G);
        return;
    end
    
    switch lower(options.WeightMode)
        case {lower('Binary')}
            if options.k > 0
                G = zeros(nSmp*(options.k+1),3);
                idNow = 0;
                for i=1:nLabel
                    classIdx = find(options.gnd==Label(i));
                    D = EuDist2(fea(classIdx,:),[],0);
                    [dump idx] = sort(D,2); % sort each row
                    clear D dump;
                    idx = idx(:,1:options.k+1);
                    
                    nSmpClass = length(classIdx)*(options.k+1);
                    G(idNow+1:nSmpClass+idNow,1) = repmat(classIdx,[options.k+1,1]);
                    G(idNow+1:nSmpClass+idNow,2) = classIdx(idx(:));
                    G(idNow+1:nSmpClass+idNow,3) = 1;
                    idNow = idNow+nSmpClass;
                    clear idx
                end
                G = sparse(G(:,1),G(:,2),G(:,3),nSmp,nSmp);
                G = max(G,G');
            else
                G = zeros(nSmp,nSmp);
                for i=1:nLabel
                    classIdx = find(options.gnd==Label(i));
                    G(classIdx,classIdx) = 1;
                end
            end
            
            if ~options.bSelfConnected
                for i=1:size(G,1)
                    G(i,i) = 0;
                end
            end
            
            W = sparse(G);
        case {lower('HeatKernel')}
            if options.k > 0
                G = zeros(nSmp*(options.k+1),3);
                idNow = 0;
                for i=1:nLabel
                    classIdx = find(options.gnd==Label(i));
                    D = EuDist2(fea(classIdx,:),[],0);
                    [dump idx] = sort(D,2); % sort each row
                    clear D;
                    idx = idx(:,1:options.k+1);
                    dump = dump(:,1:options.k+1);
                    dump = exp(-dump/(2*options.t^2));
                    
                    nSmpClass = length(classIdx)*(options.k+1);
                    G(idNow+1:nSmpClass+idNow,1) = repmat(classIdx,[options.k+1,1]);
                    G(idNow+1:nSmpClass+idNow,2) = classIdx(idx(:));
                    G(idNow+1:nSmpClass+idNow,3) = dump(:);
                    idNow = idNow+nSmpClass;
                    clear dump idx
                end
                G = sparse(G(:,1),G(:,2),G(:,3),nSmp,nSmp);
            else
                G = zeros(nSmp,nSmp);
                for i=1:nLabel
                    classIdx = find(options.gnd==Label(i));
                    D = EuDist2(fea(classIdx,:),[],0);
                    D = exp(-D/(2*options.t^2));
                    G(classIdx,classIdx) = D;
                end
            end
            
            if ~options.bSelfConnected
                for i=1:size(G,1)
                    G(i,i) = 0;
                end
            end
            
            W = sparse(max(G,G'));
        case {lower('Cosine')}
            if ~options.bNormalized
                fea = NormalizeFea(fea);
            end
            if options.k > 0
                G = zeros(nSmp*(options.k+1),3);
                idNow = 0;
                for i=1:nLabel
                    classIdx = find(options.gnd==Label(i));
                    D = fea(classIdx,:)*fea(classIdx,:)';
                    [dump idx] = sort(-D,2); % sort each row
                    clear D;
                    idx = idx(:,1:options.k+1);
                    dump = -dump(:,1:options.k+1);
                    
                    nSmpClass = length(classIdx)*(options.k+1);
                    G(idNow+1:nSmpClass+idNow,1) = repmat(classIdx,[options.k+1,1]);
                    G(idNow+1:nSmpClass+idNow,2) = classIdx(idx(:));
                    G(idNow+1:nSmpClass+idNow,3) = dump(:);
                    idNow = idNow+nSmpClass;
                    clear dump idx
                end
                G = sparse(G(:,1),G(:,2),G(:,3),nSmp,nSmp);
            else
                G = zeros(nSmp,nSmp);
                for i=1:nLabel
                    classIdx = find(options.gnd==Label(i));
                    G(classIdx,classIdx) = fea(classIdx,:)*fea(classIdx,:)';
                end
            end
            
            if ~options.bSelfConnected
                for i=1:size(G,1)
                    G(i,i) = 0;
                end
            end
            
            W = sparse(max(G,G'));
        otherwise
            error('WeightMode does not exist!');
    end
    return;
end
if bCosine && ~options.bNormalized
    Normfea = NormalizeFea(fea);
end
if strcmpi(options.NeighborMode,'KNN') && (options.k > 0)
    if ~(bCosine && options.bNormalized)
        G = zeros(nSmp*(options.k+1),3);
        for i = 1:ceil(nSmp/BlockSize)
            if i == ceil(nSmp/BlockSize)
                smpIdx = (i-1)*BlockSize+1:nSmp;
                dist = EuDist2(fea(smpIdx,:),fea,0);
                
                if bSpeed
                    nSmpNow = length(smpIdx);
                    dump = zeros(nSmpNow,options.k+1);
                    idx = dump;
                    for j = 1:options.k+1
                        [dump(:,j),idx(:,j)] = min(dist,[],2);
                        temp = (idx(:,j)-1)*nSmpNow+[1:nSmpNow]';
                        dist(temp) = 1e100;
                    end
                else
                    [dump idx] = sort(dist,2); % sort each row
                    idx = idx(:,1:options.k+1);
                    dump = dump(:,1:options.k+1);
                end
                
                if ~bBinary
                    if bCosine
                        dist = Normfea(smpIdx,:)*Normfea';
                        dist = full(dist);
                        linidx = [1:size(idx,1)]';
                        dump = dist(sub2ind(size(dist),linidx(:,ones(1,size(idx,2))),idx));
                    else
                        dump = exp(-dump/(2*options.t^2));
                    end
                end
                
                G((i-1)*BlockSize*(options.k+1)+1:nSmp*(options.k+1),1) = repmat(smpIdx',[options.k+1,1]);
                G((i-1)*BlockSize*(options.k+1)+1:nSmp*(options.k+1),2) = idx(:);
                if ~bBinary
                    G((i-1)*BlockSize*(options.k+1)+1:nSmp*(options.k+1),3) = dump(:);
                else
                    G((i-1)*BlockSize*(options.k+1)+1:nSmp*(options.k+1),3) = 1;
                end
            else
                smpIdx = (i-1)*BlockSize+1:i*BlockSize;
                dist = EuDist2(fea(smpIdx,:),fea,0);
                if bSpeed
                    nSmpNow = length(smpIdx);
                    dump = zeros(nSmpNow,options.k+1);
                    idx = dump;
                    for j = 1:options.k+1
                        [dump(:,j),idx(:,j)] = min(dist,[],2);
                        temp = (idx(:,j)-1)*nSmpNow+[1:nSmpNow]';
                        dist(temp) = 1e100;
                    end
                else
                    [dump idx] = sort(dist,2); % sort each row
                    idx = idx(:,1:options.k+1);
                    dump = dump(:,1:options.k+1);
                end
                
                if ~bBinary
                    if bCosine
                        dist = Normfea(smpIdx,:)*Normfea';
                        dist = full(dist);
                        linidx = [1:size(idx,1)]';
                        dump = dist(sub2ind(size(dist),linidx(:,ones(1,size(idx,2))),idx));
                    else
                        dump = exp(-dump/(2*options.t^2));
                    end
                end
                
                G((i-1)*BlockSize*(options.k+1)+1:i*BlockSize*(options.k+1),1) = repmat(smpIdx',[options.k+1,1]);
                G((i-1)*BlockSize*(options.k+1)+1:i*BlockSize*(options.k+1),2) = idx(:);
                if ~bBinary
                    G((i-1)*BlockSize*(options.k+1)+1:i*BlockSize*(options.k+1),3) = dump(:);
                else
                    G((i-1)*BlockSize*(options.k+1)+1:i*BlockSize*(options.k+1),3) = 1;
                end
            end
        end
        
        W = sparse(G(:,1),G(:,2),G(:,3),nSmp,nSmp);
    else
        G = zeros(nSmp*(options.k+1),3);
        for i = 1:ceil(nSmp/BlockSize)
            if i == ceil(nSmp/BlockSize)
                smpIdx = (i-1)*BlockSize+1:nSmp;
                dist = fea(smpIdx,:)*fea';
                dist = full(dist);
                if bSpeed
                    nSmpNow = length(smpIdx);
                    dump = zeros(nSmpNow,options.k+1);
                    idx = dump;
                    for j = 1:options.k+1
                        [dump(:,j),idx(:,j)] = max(dist,[],2);
                        temp = (idx(:,j)-1)*nSmpNow+[1:nSmpNow]';
                        dist(temp) = 0;
                    end
                else
                    [dump idx] = sort(-dist,2); % sort each row
                    idx = idx(:,1:options.k+1);
                    dump = -dump(:,1:options.k+1);
                end
                
                G((i-1)*BlockSize*(options.k+1)+1:nSmp*(options.k+1),1) = repmat(smpIdx',[options.k+1,1]);
                G((i-1)*BlockSize*(options.k+1)+1:nSmp*(options.k+1),2) = idx(:);
                G((i-1)*BlockSize*(options.k+1)+1:nSmp*(options.k+1),3) = dump(:);
            else
                smpIdx = (i-1)*BlockSize+1:i*BlockSize;
                dist = fea(smpIdx,:)*fea';
                dist = full(dist);
                
                if bSpeed
                    nSmpNow = length(smpIdx);
                    dump = zeros(nSmpNow,options.k+1);
                    idx = dump;
                    for j = 1:options.k+1
                        [dump(:,j),idx(:,j)] = max(dist,[],2);
                        temp = (idx(:,j)-1)*nSmpNow+[1:nSmpNow]';
                        dist(temp) = 0;
                    end
                else
                    [dump idx] = sort(-dist,2); % sort each row
                    idx = idx(:,1:options.k+1);
                    dump = -dump(:,1:options.k+1);
                end
                G((i-1)*BlockSize*(options.k+1)+1:i*BlockSize*(options.k+1),1) = repmat(smpIdx',[options.k+1,1]);
                G((i-1)*BlockSize*(options.k+1)+1:i*BlockSize*(options.k+1),2) = idx(:);
                G((i-1)*BlockSize*(options.k+1)+1:i*BlockSize*(options.k+1),3) = dump(:);
            end
        end
        W = sparse(G(:,1),G(:,2),G(:,3),nSmp,nSmp);
    end
    if bBinary
        W(logical(W)) = 1;
    end
    if isfield(options,'bSemiSupervised') && options.bSemiSupervised
        tmpgnd = options.gnd(options.semiSplit);
        
        Label = unique(tmpgnd);
        nLabel = length(Label);
        G = zeros(sum(options.semiSplit),sum(options.semiSplit));
        for idx=1:nLabel
            classIdx = tmpgnd==Label(idx);
            G(classIdx,classIdx) = 1;
        end
        Wsup = sparse(G);
        if ~isfield(options,'SameCategoryWeight')
            options.SameCategoryWeight = 1;
        end
        W(options.semiSplit,options.semiSplit) = (Wsup>0)*options.SameCategoryWeight;
    end
    if ~options.bSelfConnected
        W = W - diag(diag(W));
    end
    
    if isfield(options,'bTrueKNN') && options.bTrueKNN
    else
        W = max(W,W');
    end
    return;
end
% strcmpi(options.NeighborMode,'KNN') & (options.k == 0)
% Complete Graph
switch lower(options.WeightMode)
    case {lower('Binary')}
        error('Binary weight can not be used for complete graph!');
    case {lower('HeatKernel')}
        W = EuDist2(fea,[],0);
        W = exp(-W/(2*options.t^2));
    case {lower('Cosine')}
        W = full(Normfea*Normfea');
    otherwise
        error('WeightMode does not exist!');
end

if ~options.bSelfConnected
    for i=1:size(W,1)
        W(i,i) = 0;
    end
end

W = max(W,W');




end