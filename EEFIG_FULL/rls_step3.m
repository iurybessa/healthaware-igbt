function [K_k,theta,P]=rls_step3(P0,yk,psi_k,theta0,g,lambda)

p=size(psi_k,2);
n=size(psi_k,1);

K_k = (g*P0*psi_k')*inv(eye(n)+g*psi_k*P0*psi_k'+g*lambda*eye(n));
theta = theta0 + K_k*(yk-(psi_k*theta0));
P = (1/lambda)*(P0-P0*psi_k'*inv(psi_k*P0*psi_k'+lambda*eye(n))*psi_k*P0);
            if any(isnan(theta))
                P
            end
end
