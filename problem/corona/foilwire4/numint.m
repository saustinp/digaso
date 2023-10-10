function [SaPeek,SPeek] = numint(e, eta, e0, ind)

M = size(e,2);
SPeek = 0;
SaPeek = zeros(M,1);
for m = 1:M        
    ind0 = ind(m,:);
    ind0 = ind0(ind0>0);
    if isempty(ind0)==0
    if length(ind0)==1
        if ind0>1        
            SaPeek(m) = 0.5*sum((e(2:ind0,m)+e(1:(ind0-1),m)-2*e0).*(eta(2:ind0,m)-eta(1:(ind0-1),m)));
        end        
    else
        N = length(ind0);
        SPeek = zeros(N,1);
        for n = 1:N
            if n == 1
                SPeek(n) = 0.5*sum((e(2:ind0(n),m)+e(1:(ind0(n)-1),m)-2*e0).*(eta(2:ind0(n),m)-eta(1:(ind0(n)-1),m)));
            else
                SPeek(n) = 0.5*sum((e((ind0(n-1)+1):ind0(n),m)+e(ind0(n-1):(ind0(n)-1),m)-2*e0).*(eta((ind0(n-1)+1):ind0(n),m)-eta(ind0(n-1):(ind0(n)-1),m)));
            end
        end    
        in = SPeek>=0;
        SaPeek(m) = sum(SPeek(in));        
    end
    end
end


% function SaPeek = numint(e, eta, e0, eta0, ind0)
% 
% if isempty(eta0)
%     SaPeek = 0;
% elseif length(eta0)==1
% %     SaPeek = 0;
% %     for i = 1:(ind0-1)
% %         SaPeek = SaPeek + 0.5*(e(i)+e(i+1)-2*e0)*(eta(i+1)-eta(i));
% %     end
%     SaPeek = 0.5*sum((e(2:ind0)+e(1:(ind0-1))-2*e0).*(eta(2:ind0)-eta(1:(ind0-1))));
% else   
%     N = length(eta0);
%     SPeek = zeros(N,1);
%     for n = 1:N
%         if n == 1
%             SPeek(n) = 0.5*sum((e(2:ind0(n))+e(1:(ind0(n)-1))-2*e0).*(eta(2:ind0(n))-eta(1:(ind0(n)-1))));
%         else
%             SPeek(n) = 0.5*sum((e((ind0(n-1)+1):ind0(n))+e(ind0(n-1):(ind0(n)-1))-2*e0).*(eta((ind0(n-1)+1):ind0(n))-eta(ind0(n-1):(ind0(n)-1))));
%         end
%     end    
%     ind = SPeek>=0;
%     SaPeek = sum(SPeek(ind));
% end
% 
% 
