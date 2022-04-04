function out=mapek(r,rhat,k,iEOL)

MAX=100;
MIN=0;
r2=min(MAX,max(MIN,rhat));

mapesum=0;
for i=k:iEOL
    mapesum=mapesum+abs(r(i)-r2(i))/r(i);
end


out=(100/(iEOL-k+1))*mapesum;
% if isnan(out)
%     out=(100/(k-N+1))*sum(abs(r2(k+N:k+N)-rhat(k+1:k+N))./r(k+1:k+N));
% end
end