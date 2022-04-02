function out=rak(rhat,r,k)

MAX=100;
MIN=0;
r2=min(MAX,max(MIN,rhat));
out=max(1-abs(r(k)-r2(k))./r(k),0);
% if isnan(out)
%     out=(100/(k-N+1))*sum(abs(r2(k+N:k+N)-rhat(k+1:k+N))./r(k+1:k+N));
% end
end