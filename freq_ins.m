function f_ins = freq_ins(X, f)

% X: time-frequency function
% f: f-axis tick
%repmat(堆疊矩陣)
[~,T]=size(X);
%縱向排列f,排T長度
f_ins=sum(repmat(f(:),1,T).*(abs(X)))./(sum(abs(X))); %小心分母為零
f_ins(isnan(f_ins))=0; %isnan為在無限大的地方呈現1 ; 使其變為零

end