function ind = bisection(e, e0, sop)

[N,M] = size(e);
ind = zeros(M,1);

tmp = zeros(1,1);
if sop==1
    for m = 1:M        
        k = 0;
        for n = (N-1):-1:1
            if (e(n,m)>=e0) && (e(n+1,m)<=e0)
                k = k + 1;
                if abs(e(n,m)-e0)<abs(e(n+1,m)-e0)                
                    tmp(k) = n;
                else                
                    tmp(k) = n+1; 
                end
                %break;                
            end        
        end        
        tmq = unique(tmp(1:k));
        ind(m,1:length(tmq)) = tmq;
    end
else
    for m = 1:M
        k = 0;
        for n = 1:1:N-1            
            if (e(n,m)>=e0) && (e(n+1,m)<=e0)
                k = k + 1;
                if abs(e(n,m)-e0)<abs(e(n+1,m)-e0)                
                    tmp(k) = n;
                else                
                    tmp(k) = n+1; 
                end
                %break;                
            end        
        end        
        tmq = unique(tmp(1:k));
        ind(m,1:length(tmq)) = tmq;        
    end    
end

return;


% function [eta0,i0] = bisection(e, eta, e0)
% 
% emin = min(e);
% emax = max(e);
% 
% i0 = [];
% if (emin<e0) && (e0<emax)
%     N = length(e);    
%     for i = 1:N-1
%         if (e(i)>=e0) && (e(i+1)<=e0)
%             if abs(e(i)-e0)<abs(e(i+1)-e0)                
%                 i0 = [i0 i];
%             else                
%                 i0 = [i0 i+1];
%             end
%         end
%     end
% end
% i0 = unique(i0);
% eta0 = eta(i0);
% 
% 
% return;
% 
