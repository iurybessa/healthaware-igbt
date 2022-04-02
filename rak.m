function out=rak(r,rhat,k)

MAX=100;
MIN=0;
r2=min(MAX,max(MIN,r));
out=1-abs(r2(k)-rhat(k))./r(k);
% if isnan(out)
%     out=(100/(k-N+1))*sum(abs(r2(k+N:k+N)-rhat(k+1:k+N))./r(k+1:k+N));
% end
end