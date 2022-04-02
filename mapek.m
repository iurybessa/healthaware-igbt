function out=mapek(r,rhat,k)

MAX=100;
MIN=0;
r2=min(MAX,max(MIN,r));
N=min(find(r2<=0));
if(isempty(N))
    N=1;
end
out=(100/(N))*sum(abs(r2(k+1:k+N)-rhat(k+1:k+N))./r(k+1:k+N));
% if isnan(out)
%     out=(100/(k-N+1))*sum(abs(r2(k+N:k+N)-rhat(k+1:k+N))./r(k+1:k+N));
% end
end