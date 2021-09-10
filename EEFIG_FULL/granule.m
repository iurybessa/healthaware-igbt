classdef granule
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        counter
        C
        m
        esm
        is_anomaly
        ignore
        anomaly_points
        beta
        alpha
        L_it
        R_it
        m_it
        a_it
        b_it
        a
        b
        L
        R
        k
        w
        idx
        Q
        Q_it
        gsum
        a1_it
        b1_it
        L1_it
        R1_it
        m1_it
        is_improved
        A
        B
        lambda
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        %% Granule Class
        function [obj] = granule(F)
            if nargin ~=0
                m = F(2);
                obj(m) = granule;
                for i = 1:m
                    obj(i) = obj(i).gran_init(F(1));
                end
            end
        end % end_granule
        %% Initial Granule Evolution
        function [gran,is_anomaly] = evolv_gran_init(prototype,g,xk,t1,varargin)
            global r
            if (sum(isnan(xk))>0)
                gran = prototype;
                is_anomaly = 0;
                return;
            end
            mahal_xk = (xk-prototype.m)*(prototype.C)*(xk-prototype.m)';
            mahal1_xk = (xk-prototype.m)*(prototype.C)*(xk-prototype.m)';
            is_anomaly = 0;
            
            gran = prototype;
            gran.beta = prototype.beta + (g^2);
            gran.alpha = prototype.alpha + g;
            
            gama = prototype.alpha * (((gran.alpha)^2) - (gran.beta));
            gama = (gama) / ((gran.alpha) *(((prototype.alpha)^2) - gran.beta));
            
            ForceUpdate = 1;
            if (nargin>3)
                ForceUpdate = varargin{1};
            end
            
            if (mahal_xk > t1)
                
                if (mahal_xk>chi2inv(0.999,size(prototype.m,2))) && (ForceUpdate>1)
                    ForceUpdate = 0;
                end
                
                is_anomaly = 1;
                
                if ForceUpdate == 0
                    return;
                end
                
            end
            
            delta = (gran.alpha) * (((prototype.alpha)^2) - prototype.beta);
            delta = (delta) / (((prototype.alpha) * g) * (g+(gran.alpha)-2));
            mahal_xk = (mahal_xk) + (delta);
            gran.m = (prototype.m) + ((g/gran.gsum) * (xk-prototype.m)) ;
            
            if (sum(isnan(gran.m))>0)
                gran = prototype;
                is_anomaly = 0;
                return;
            end
            
            gran.C = (prototype.C) - (((gran.C * (xk-prototype.m)') * (xk-prototype.m) * (prototype.C)) / (mahal_xk));
            gran.C = (gama) * (gran.C);
            
            gran.counter = gran.counter + 1;
            
            idx1 = gran.idx + 1;
            
            for ii = 1:size(xk,1)
                for jj = 1:size(xk,2)
                    gran.m1_it(ii,jj) = ((idx1-1)/(idx1))*(gran.m(ii,jj))+(1/idx1)*xk(ii,jj);
                    check_L = 0;
                    check_R = 0;
                    if xk(ii,jj) < gran.m1_it(ii,jj)
                        kk1(ii,jj) = gran.k(ii,jj)+1;
                        ww1(ii,jj) = gran.w(ii,jj);
                        gran.L1_it(ii,jj) = (gran.L(ii,jj)*(kk1(ii,jj)-1) + (gran.m1_it(ii,jj)-xk(ii,jj)))/kk1(ii,jj);
                        check_L = 1;
                    else
                        kk1(ii,jj) = gran.k(ii,jj);
                        ww1(ii,jj) = gran.w(ii,jj)+1;
                        gran.R1_it(ii,jj) = (gran.R(ii,jj)*(ww1(ii,jj)-1) + (-gran.m1_it(ii,jj)+xk(ii,jj))) / ww1(ii,jj);
                        check_R = 1;
                    end
                    if check_L == 0
                        gran.L1_it(ii,jj) = gran.L(ii,jj);
                    end
                    if check_R == 0
                        gran.R1_it(ii,jj) = gran.R(ii,jj);
                    end
                    gran.a1_it(ii,jj) = max(gran.m1_it(ii,jj)-4.*gran.L1_it(ii,jj),0)+(gran.m1_it(ii,jj)-max(gran.m1_it(ii,jj)-4.*gran.L1_it(ii,jj),0))*0;
                    gran.b1_it(ii,jj) = min(gran.m1_it(ii,jj)+4.*gran.R1_it(ii,jj),Inf)-(min(gran.m1_it(ii,jj)+4.*gran.R1_it(ii,jj),Inf)-gran.m1_it(ii,jj))*0;
                end
            end
            gran.L_it = [gran.L_it ; gran.L1_it];
            gran.L = gran.L_it(end,:);
            gran.R_it = [gran.R_it ; gran.R1_it];
            gran.R = gran.R_it(end,:);
            gran.a_it = [gran.a_it ; gran.a1_it];
            gran.a = gran.a_it(end,:);
            gran.b_it = [gran.b_it ; gran.b1_it];
            gran.b = gran.b_it(end,:);
            gran.k = kk1;
            gran.w = ww1;
            gran.idx = idx1;
            gran.m_it = [gran.m_it ; gran.m1_it];
            gran.m  = gran.m_it(end,:);
            gran.gsum = gran.gsum+g;
            Q_xk = (mahal1_xk) * (gran.gsum); % h.c
            %             Q_xk = (mahal1_xk) * g; % h.c
            gran.Q_it = [gran.Q_it ; Q_xk];
            gran.Q = cumsum(gran.Q_it);
        end %end_evolv_gran_init
        %% Granule initialization
        function [gran] = gran_init(oldrule,p,varargin)
            global r
            gran = granule;
            if (nargin>2)
                gran = oldrule;
                gran.m_it = mean(varargin{1}(1:2,1:p));
                for j = 1:size(varargin{1},2)
                    gran.L_it(1,j) = min(varargin{1}(:,j));
                    gran.R_it(1,j) = max(varargin{1}(:,j));
                end
                for j = 1:size(varargin{1},2)
                    gran.a_it(1,j) = max(gran.m_it(1,j)-4.*gran.L_it(1,j),0)+(gran.m_it(1,j)-max(gran.m_it(1,j)-4.*gran.L_it(1,j),0))*0;
                    gran.b_it(1,j) = min(gran.m_it(1,j)+4.*gran.R_it(1,j),Inf)-(min(gran.m_it(1,j)+4.*gran.R_it(1,j),Inf)-gran.m_it(1,j))*0;
                end
                kk = ones(1,size(varargin{1},2));
                ww = ones(1,size(varargin{1},2));
                gran.a = gran.a_it(end,:);
                gran.b = gran.b_it(end,:);
                gran.L = gran.L_it(end,:);
                gran.R = gran.R_it(end,:);
                gran.k = kk(end,:);
                gran.w = ww(end,:);
                gran.m = gran.m_it(end,:);
                gran.idx = 2;
                gran.counter = 2;
                gran.alpha = 2;
                gran.beta = 2;
                gran.gsum = 2;
                Q_xk1 = ((varargin{1}(1,1:p)-gran.m)*gran.C*(varargin{1}(1,1:p)-gran.m)')*1;
                Q_xk2 = ((varargin{1}(2,1:p)-gran.m)*gran.C*(varargin{1}(2,1:p)-gran.m)')*gran.gsum;
                gran.Q_it = [Q_xk1;Q_xk2];
                gran.Q = cumsum(gran.Q_it);
                gran.A = [];
                gran.B = [];
                
                for i = 3:1:size(varargin{1},1)
                    t1 = chi2inv(0.999,p);
                    [g,~,~,~] = data_evaluation(gran,varargin{1}(i,1:p),t1);
                    [gran,~] = gran.evolv_gran_init(g,varargin{1}(i,1:p),t1,1);
                end
            else
                gran.counter = 0;
                gran.alpha = 0;
                gran.beta = 0;
                gran.C = eye(p);
                gran.m = zeros(1,p);
                gran.esm = 0;
                gran.is_anomaly = 0;
                gran.ignore = 0;
                gran.anomaly_points = 0;
                gran.L_it = zeros(1,p);
                gran.R_it = zeros(1,p);
                gran.m_it = zeros(1,p);
                gran.a_it = zeros(1,p);
                gran.b_it = zeros(1,p);
                gran.a = zeros(1,p);
                gran.b = zeros(1,p);
                gran.L = zeros(1,p);
                gran.R = zeros(1,p);
                gran.k = zeros(1,p);
                gran.w = zeros(1,p);
                gran.idx = 0;
                gran.alpha = 0;
                gran.beta = 0;
                gran.Q_it = 0;
                gran.Q = 0;
            end
        end %end_gran_init
        %% Granule Saving
        function [gran] = Duplicate (oldrule)
            gran = oldrule.gran_init(numel(oldrule.m));
            gran.counter = oldrule.counter;
            gran.alpha = 0;
            gran.beta = 0;
            gran.m = oldrule.m;
            gran.C = oldrule.C;
            gran.L_it = oldrule.L_it;
            gran.R_it = oldrule.R_it;
            gran.m_it = oldrule.m_it;
            gran.a_it = oldrule.a_it;
            gran.b_it = oldrule.b_it;
            gran.L1_it = oldrule.L1_it;
            gran.R1_it = oldrule.R1_it;
            gran.m1_it = oldrule.m1_it;
            gran.a1_it = oldrule.a1_it;
            gran.b1_it = oldrule.b1_it;
            gran.a = oldrule.a;
            gran.b = oldrule.b;
            gran.L = oldrule.L;
            gran.R = oldrule.R;
            gran.k = oldrule.k;
            gran.w = oldrule.w;
            gran.idx = oldrule.idx;
            gran.Q_it = oldrule.Q_it;
            gran.Q = oldrule.Q;
        end % end duplicate
        %% Granule Evolution
        function [gran,is_anomaly] =  evolv_gran(prototype,xk,t1,g,wk_i,wk,...
                l0_i,l0,l1_i,l1,f_i,f,varargin)
            global r
            if (sum(isnan(xk))>0)
                gran = prototype;
                is_anomaly = 0;
                return;
            end
            mahal_xk = (xk-prototype.m)*(prototype.C)*(xk-prototype.m)';
            mahal1_xk = (xk-prototype.m)*(prototype.C)*(xk-prototype.m)';
            is_anomaly = 0;
            
            gran = prototype;
            gran_test = gran;
            gran.beta = prototype.beta + (g^2);
            gran.alpha = prototype.alpha + g;
            
            gama = prototype.alpha * (((gran.alpha)^2) - (gran.beta));
            gama = gama / ((gran.alpha) *(((prototype.alpha)^2) - gran.beta));
            
            ForceUpdate = 1;
            if (nargin>3)
                ForceUpdate = varargin{1};
            end
            
            if (mahal_xk > t1)
                
                if (mahal_xk>chi2inv(0.999,size(prototype.m,2))) && (ForceUpdate>1)
                    ForceUpdate = 0;
                end
                
                is_anomaly = 1;
                
                if ForceUpdate == 0
                    return;
                end
                
            end
            
            delta = (gran.alpha) * (((prototype.alpha)^2) - prototype.beta);
            delta = (delta) / (((prototype.alpha) * g) * (g+(gran.alpha)-2));
            mahal_xk = (mahal_xk) + (delta);
            gran.m = (prototype.m) + ((g/gran.gsum) * (xk-prototype.m)) ;
            %             gran.m = prototype.m+(g)*(xk-prototype.m);
            if (sum(isnan(gran.m))>0)
                gran = prototype;
                is_anomaly = 0;
                return;
            end
            
            gran.C = (prototype.C) - (((gran.C * (xk-prototype.m)') * (xk-prototype.m) * (prototype.C)) / (mahal_xk));
            gran.C = (gama) * (gran.C);
            
            kk1 = gran.k;
            ww1 = gran.w;
            
            idx1 = gran.idx+1;
            
            for ii = 1:size(xk,1)
                for jj = 1:size(xk,2)
                    gran.m1_it(ii,jj) = ((idx1-1)/(idx1))*(gran.m(ii,jj))+(1/idx1)*xk(ii,jj);
                    check_L = 0;
                    check_R = 0;
                    if xk(ii,jj) < gran.m1_it(ii,jj)
                        kk1(ii,jj) = gran.k(ii,jj)+1;
                        ww1(ii,jj) = gran.w(ii,jj);
                        gran.L1_it(ii,jj) = (gran.L(ii,jj)*(kk1(ii,jj)-1) + (gran.m1_it(ii,jj)-xk(ii,jj)))/kk1(ii,jj);
                        check_L = 1;
                    else
                        kk1(ii,jj) = gran.k(ii,jj);
                        ww1(ii,jj) = gran.w(ii,jj)+1;
                        gran.R1_it(ii,jj) = (gran.R(ii,jj)*(ww1(ii,jj)-1) + (-gran.m1_it(ii,jj)+xk(ii,jj))) / ww1(ii,jj);
                        check_R = 1;
                    end
                    if check_L == 0
                        gran.L1_it(ii,jj) = gran.L(ii,jj);
                    end
                    if check_R == 0
                        gran.R1_it(ii,jj) = gran.R(ii,jj);
                    end
                    gran.a1_it(ii,jj) = max(gran.m1_it(ii,jj)-4.*gran.L1_it(ii,jj),0)+(gran.m1_it(ii,jj)-max(gran.m1_it(ii,jj)-4.*gran.L1_it(ii,jj),0))*0;
                    gran.b1_it(ii,jj) = min(gran.m1_it(ii,jj)+4.*gran.R1_it(ii,jj),Inf)-(min(gran.m1_it(ii,jj)+4.*gran.R1_it(ii,jj),Inf)-gran.m1_it(ii,jj))*0;
                end
            end
            gran.L_it = [gran.L_it ; gran.L1_it];
            gran.L = gran.L_it(end,:);
            gran.R_it = [gran.R_it ; gran.R1_it];
            gran.R = gran.R_it(end,:);
            gran.a_it = [gran.a_it ; gran.a1_it];
            gran.a = gran.a_it(end,:);
            gran.b_it = [gran.b_it ; gran.b1_it];
            gran.b = gran.b_it(end,:);
            gran.k = kk1;
            gran.w = ww1;
            gran.idx = idx1;
            gran.m_it = [gran.m_it ; gran.m1_it];
            gran.m  = gran.m_it(end,:);
            [gran,is_Q_improved] = granularity_allocation(gran,gran_test,...
                mahal1_xk,g,wk_i,wk,l0_i,l0,l1_i,l1,f_i,f);
            
            if is_Q_improved == 1
                gran.counter = gran.counter + 1;
                gran.is_improved = is_Q_improved;
            elseif is_Q_improved ==0
                gran.counter = gran.counter;
                gran.is_improved = is_Q_improved;
            end
            
        end %end_evolv_gran
        %% Granule Step
        function [newprototype,is_anomaly] = gran_step(prototype,xk,sp,t1,wk_i,wk,...
                lambda0_i,lambda0,lambda1_i,lambda1,f_i,f,varargin)
            newprototype = prototype;
            g = 1;
            if (nargin>3)
                g = varargin{1};
            end
            if (prototype.counter<sp)
                [newprototype,is_anomaly] = newprototype.evolv_gran(xk,t1,g,wk_i,wk,...
                    lambda0_i,lambda0,lambda1_i,lambda1,f_i,f,1);
            else
                [newprototype,is_anomaly] = newprototype.evolv_gran(xk,t1,g,wk_i,wk,...
                    lambda0_i,lambda0,lambda1_i,lambda1,f_i,f,2);
            end
        end %end_granule_step
    end %end_methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %end_classdef

