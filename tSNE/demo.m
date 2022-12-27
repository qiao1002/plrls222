clc 
clear
%% load data sets ('in_X' and 'true_labs')
dataset = {'Darmanis','Kold','Tasic','Zeisel','Ramskold',	'islet',	'Treutlein','Ting','Goolam','Deng','19_Engel','Pollen'}
% dataset = {'Ginhoux','Grover','Leng','Buettner','mECS','Zheng','Human','Bre','Zeisel','Ramskold','Deng','Darmanis','Tasic','Usoskin','islet','Treutlein','Ting','Kold','19_Engel','Pollen','Goolam', 'Kold','Grover',}; %% four datasets tested on the paper
for num=4:1:4
    load(['Data_',dataset{num},'.mat']);
    %% 预处理
%      [in_X,] = FilterGenesZero(in_X);
     fea=double(in_X);
     gnd=true_labs(:);
     num_class =length(unique(gnd)); % The number of classes 
      %% spectral clustering
% %      SimGraph= func_simlarity(double(in_X));%高斯核
%      similarity=corr(double(in_X)',double(in_X)','type','pearson');%皮尔逊
%     [result_label, kerNS]= SpectralClustering(similarity,num_class);  %%
%     % ---------- Evaluation measurement  -------- %
%     NMI=Cal_NMI_newused(gnd, result_label);
%     Purity = purity(max(gnd),result_label, gnd);
%     ARI=Contingency_ARI_newused(gnd, result_label);
     for  jjj=1:30
        mappedX = tsne(double(in_X), true_labs, num_class );
        ind_tsne = cal_new_clus(mappedX,num_class );
        nmi_tsne=Cal_NMI(ind_tsne,true_labs);
        ari_tsne=Cal_ARI(ind_tsne,true_labs);
        [result1(jjj,:)]=nmi_tsne;
        [result2(jjj,:)]=ari_tsne;
%         [result3(jjj,:)]=Purity;
    end
%    [result1(num,:)]=NMI;
%     [result2(num,:)]=ARI;
     [result1_avg(num,:)]=mean(result1);
    [result2_avg(num,:)]=mean(result2);
end
% nmi_average=mean(result1)
% ari_average=mean(result2)



