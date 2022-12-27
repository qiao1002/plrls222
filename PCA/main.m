clear all
close all
dataset = {'Darmanis','Kold','Tasic','Zeisel','Ramskold',	'islet',	'Treutlein','Ting','Goolam','Deng','19_Engel','Pollen'};
for num=1:1:12
    load(['Data_',dataset{num},'.mat']);
 
    [in_X,] = FilterGenesZero(in_X);
    X=in_X';
    maxIter = 50;
    [m,n]=size(X);%%%
    beta =0.56;
    k=numel(unique(true_labs));
    fea = X;
    options = [];
    options.WeightMode = 'Binary';
    W = constructW(fea',options); 
    [U,V,err]= PCA(X,W,beta,k);
    cluster_result();
    [result1(num,:)]=nmi_mn;
    [result2(num,:)]=ari_mn;
end
%=============  evaluate AC: accuracy ==============
% label2 = kmeans(V,nClass);
% res = bestMap(gnd,label2);
% AC21=length(find(gnd == res))/length(gnd);
% MIhat2 = MutualInfo(gnd,label2);%判断聚类效果
% disp(['Clustering in the NMF AC11 and MI: ',num2str(AC21),', ',num2str(MIhat2)]);
%end
%=======================特征选择============
% u = sum(abs(U(:,1:k)),2);
% S=sort(u,'descend');
% number = 500;
% [index, num] =  selectfeature(u, number);
% filename = ['PAAD_CHOL_COAD result 500(50次).txt'];
% fid = fopen(filename, 'wt');
% 
% fprintf(fid,'number = %d\n', num);
% for i = 1 : num
%     fprintf(fid,'%s\n', PAAD_CHOL_COAD_featurename3{index(i)});
% end
% fclose(fid);





