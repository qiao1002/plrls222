function [score]=fraction(A,sel,yn)

[m,n]=size(A);
if sel==0
    D =pdist(A,'minkowski',yn);%yn=1-Âü¹þ¶Ù£»yn=2-Å·Ê½¾àÀë
else
    D =pdist(A,'chebychev');%ÇÐ±ÈÑ©·ò
end
%D =pdist(A,'correlation');
% score = squareform(D); 
score=zeros(n,n);
m=1;
for i=1:n
    for j=i+1:n
        score(i,j)=D(1,m);
        m=m+1;
        i;
    end
end
score=score./max(max(score));
score=score+score';
score=score./max(max(score));
score=1./(score+1);
% squareform(D)

end

% 	num=15;
% 	[U,S,V] = svd(X);
% 	l = min((m,n))
% 	sigma = zeros(l,l);
% 	for i=0:1:num
% 		sigma[i,i] = S[i,i]
%     end
%     new_expression=U*sigma*V;
% 	return(new_expression)
%     
%     
%   %%Ï¡ÊèÐÔ
%     [m,n] = size(X);
% 	zeros_num = sum(X == 0)
% 	sparsity = zeros_num/(n*m)
% 	return(sparsity)
% 
