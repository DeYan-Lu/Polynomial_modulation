function f_ins = freq_ins(X, f)

% X: time-frequency function
% f: f-axis tick
%repmat(���|�x�})
[~,T]=size(X);
%�a�V�ƦCf,��T����
f_ins=sum(repmat(f(:),1,T).*(abs(X)))./(sum(abs(X))); %�p�ߤ������s
f_ins(isnan(f_ins))=0; %isnan���b�L���j���a��e�{1 ; �Ϩ��ܬ��s

end