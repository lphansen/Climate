function res = finiteDifference(value, dimension, order, dlt, scheme, capbound)
% 
% if sum(size(dlt))>2
%     s = size(dlt);
%     if s(1,2)>1
%         dlt=dlt';
%     end
% end
    
res = zeros(size(value));
if nargin > 5
    cap = capbound;
else
    cap = nan;
end

if sum(size(dlt))<=2 % if step size dlt is fixed
    if strcmp(scheme, 'central')
        if order == 1   % first order derivatives
            if dimension == 1
                res(2:end-1,:,:) = (1 ./ (2.*dlt)) .* ( value(3:end,:,:) - value(1:end-2,:,:));
                res(end,:,:) = (1 ./ dlt) .* (value(end,:,:) - value(end-1,:,:));
                res(1,:,:) = (1 ./ dlt) .* (value(2,:,:) - value(1,:,:));
            elseif dimension == 2
                res(:,2:end-1,:) = (1 ./ (2.*dlt)) .* ( value(:,3:end,:) - value(:,1:end-2,:));
                res(:,end,:) = (1 ./ dlt) .* (value(:,end,:) - value(:,end-1,:));
                res(:,1,:) = (1 ./ dlt) .* (value(:,2,:) - value(:,1,:));
            else
                res(:,:,2:end-1) = (1 ./ (2.*dlt)) .* ( value(:,:,3:end) - value(:,:,1:end-2));
                res(:,:,end) = (1 ./ dlt) .* (value(:,:,end) - value(:,:,end-1));
                res(:,:,1) = (1 ./ dlt) .* (value(:,:,2) - value(:,:,1));
            end
            
        elseif order == 2
            if dimension == 1
                res(2:end-1,:,:) = (1 ./ (dlt)^2) .* ( value(3:end,:,:) + value(1:end-2,:,:) - 2 .* value(2:end-1,:,:));
                res(end,:,:) = (1 ./ dlt ^ 2) .* (value(end,:,:) + value(end-2,:,:) - 2 .* value(end-1,:,:));
                res(1,:,:) = (1 ./ dlt^2) .* (value(3,:,:) + value(1,:,:) - 2 .* value(2,:,:));
            elseif dimension == 2
                res(:,2:end-1,:) = (1 ./ (dlt)^2) .* ( value(:,3:end,:) + value(:,1:end-2,:) - 2 .* value(:,2:end-1,:));
                res(:,end,:) = (1 ./ dlt ^ 2) .* (value(:,end,:) + value(:,end-2,:) - 2 .* value(:,end-1,:));
                res(:,1,:) = (1 ./ dlt^2) .* (value(:,3,:) + value(:,1,:) - 2 .* value(:,2,:));
            else
                res(:,:,2:end-1) = (1 ./ (dlt)^2) .* ( value(:,:,3:end) + value(:,:,1:end-2) - 2 .* value(:,:,2:end-1));
                res(:,:,end) = (1 ./ dlt ^ 2) .* (value(:,:,end) + value(:,:,end-2) - 2 .* value(:,:,end-1));
                res(:,:,1) = (1 ./ dlt^2) .* (value(:,:,3) + value(:,:,1) - 2 .* value(:,:,2));
            end
        else
            disp('order not supported');
        end
        
    else
        if order == 1   % first order derivatives
            if dimension == 1
                res(2:end-1,:,:) = (1 ./ (dlt)) .* ( value(3:end,:,:) - value(2:end-1,:,:));
                res(end,:,:) = (1 ./ dlt) .* (value(end,:,:) - value(end-1,:,:));
                res(1,:,:) = (1 ./ dlt) .* (value(2,:,:) - value(1,:,:));
            elseif dimension == 2
                res(:,2:end-1,:) = (1 ./ (dlt)) .* ( value(:,3:end,:) - value(:,2:end-1,:));
                res(:,end,:) = (1 ./ dlt) .* (value(:,end,:) - value(:,end-1,:));
                res(:,1,:) = (1 ./ dlt) .* (value(:,2,:) - value(:,1,:));
            else
                res(:,:,2:end-1) = (1 ./ (dlt)) .* ( value(:,:,3:end) - value(:,:,2:end-1));
                res(:,:,end) = (1 ./ dlt) .* (value(:,:,end) - value(:,:,end-1));
                res(:,:,1) = (1 ./ dlt) .* (value(:,:,2) - value(:,:,1));
            end
            
        elseif order == 2
            if dimension == 1
                res(2:end-1,:,:) = (1 ./ (dlt)^2) .* ( value(3:end,:,:) + value(1:end-2,:,:) - 2 .* value(2:end-1,:,:));
                res(end,:,:) = (1 ./ dlt ^ 2) .* (value(end,:,:) + value(end-2,:,:) - 2 .* value(end-1,:,:));
                res(1,:,:) = (1 ./ dlt^2) .* (value(3,:,:) + value(1,:,:) - 2 .* value(2,:,:));
            elseif dimension == 2
                res(:,2:end-1,:) = (1 ./ (dlt)^2) .* ( value(:,3:end,:) + value(:,1:end-2,:) - 2 .* value(:,2:end-1,:));
                res(:,end,:) = (1 ./ dlt ^ 2) .* (value(:,end,:) + value(:,end-2,:) - 2 .* value(:,end-1,:));
                res(:,1,:) = (1 ./ dlt^2) .* (value(:,3,:) + value(:,1,:) - 2 .* value(:,2,:));
            else
                res(:,:,2:end-1) = (1 ./ (dlt)^2) .* ( value(:,:,3:end) + value(:,:,1:end-2) - 2 .* value(:,:,2:end-1));
                res(:,:,end) = (1 ./ dlt ^ 2) .* (value(:,:,end) + value(:,:,end-2) - 2 .* value(:,:,end-1));
                res(:,:,1) = (1 ./ dlt^2) .* (value(:,:,3) + value(:,:,1) - 2 .* value(:,:,2));
            end
        else
            disp('order not supported');
        end        
        
    end
        
else  % if step size dlt is a vector
    if strcmp(scheme, 'central')
        if order == 1   % first order derivatives
            if dimension == 1
                dlt_multi = reshape(dlt,[length(dlt),1,1]);
                res(2:end-1,:,:) = (1 ./ (dlt_multi(2:end) + dlt_multi(1:end-1))) .* ( value(3:end,:,:) - value(1:end-2,:,:));
                res(end,:,:) = (1 ./ dlt(end)) .* (value(end,:,:) - value(end-1,:,:));
                res(1,:,:) = (1 ./ dlt(1)) .* (value(2,:,:) - value(1,:,:));
            elseif dimension == 2
                dlt_multi = reshape(dlt,[1,length(dlt),1]);
                res(:,2:end-1,:) = (1 ./ (dlt_multi(2:end) + dlt_multi(1:end-1))) .* ( value(:,3:end,:) - value(:,1:end-2,:));
                res(:,end,:) = (1 ./ dlt(end)) .* (value(:,end,:) - value(:,end-1,:));
                res(:,1,:) = (1 ./ dlt(1)) .* (value(:,2,:) - value(:,1,:));
            else
                dlt_multi = reshape(dlt,[1,1,length(dlt)]);
                res(:,:,2:end-1) = (1 ./ (dlt_multi(2:end) + dlt_multi(1:end-1))) .* ( value(:,:,3:end) - value(:,:,1:end-2));
                res(:,:,end) = (1 ./ dlt(end)) .* (value(:,:,end) - value(:,:,end-1));
                res(:,:,1) = (1 ./ dlt(1)) .* (value(:,:,2) - value(:,:,1));
            end
            
        elseif order == 2
            if dimension == 1
                dlt_multi = reshape(dlt,[length(dlt),1,1]);
                res(2:end-1,:,:) = ((1 ./ dlt_multi(2:end)) .* ( value(3:end,:,:) - value(2:end-1,:,:)) - (1 ./ dlt_multi(1:end-1)) .* ( value(2:end-1,:,:) - value(1:end-2,:,:))) ./ dlt_multi(1:end-1);
                res(end,:,:) = (1 ./ dlt(end) ^ 2) .* (value(end,:,:) + value(end-2,:,:) - 2 .* value(end-1,:,:));
                res(1,:,:) = (1 ./ dlt(1)^2) .* (value(3,:,:) + value(1,:,:) - 2 .* value(2,:,:));
            elseif dimension == 2
                dlt_multi = reshape(dlt,[1,length(dlt),1]);
                res(:,2:end-1,:) = (1 ./ dlt_multi(2:end) .* ( value(:,3:end,:) - value(:,2:end-1,:)) - 1 ./ dlt_multi(1:end-1) .* ( value(:,2:end-1,:) - value(:,1:end-2,:))) ./ dlt_multi(1:end-1);
                res(:,end,:) = (1 ./ dlt(end) ^ 2) .* (value(:,end,:) + value(:,end-2,:) - 2 .* value(:,end-1,:));
                res(:,1,:) = (1 ./ dlt(1)^2) .* (value(:,3,:) + value(:,1,:) - 2 .* value(:,2,:));
            else
                dlt_multi = reshape(dlt,[1,1,length(dlt)]);
                res(:,:,2:end-1) = (1 ./ dlt_multi(2:end) .* ( value(:,:,3:end) - value(:,:,2:end-1)) - 1 ./ dlt_multi(1:end-1) .* ( value(:,:,2:end-1) - value(:,:,1:end-2))) ./ dlt_multi(1:end-1);
                res(:,:,end) = (1 ./ dlt(end) ^ 2) .* (value(:,:,end) + value(:,:,end-2) - 2 .* value(:,:,end-1));
                res(:,:,1) = (1 ./ dlt(1)^2) .* (value(:,:,3) + value(:,:,1) - 2 .* value(:,:,2));
            end
        else
            disp('order not supported');
        end
        
    else
        if order == 1   % first order derivatives
            if dimension == 1
                res(2:end-1,:,:) = 1 ./ dlt(2:end) .* ( value(3:end,:,:) - value(1:end-2,:,:));
                res(end,:,:) = (1 ./ dlt(end)) .* (value(end,:,:) - value(end-1,:,:));
                res(1,:,:) = (1 ./ dlt(1)) .* (value(2,:,:) - value(1,:,:));
            elseif dimension == 2
                res(:,2:end-1,:) = 1 ./ dlt(2:end) .* ( value(:,3:end,:) - value(:,1:end-2,:));
                res(:,end,:) = (1 ./ dlt(end)) .* (value(:,end,:) - value(:,end-1,:));
                res(:,1,:) = (1 ./ dlt(1)) .* (value(:,2,:) - value(:,1,:));
            else
                res(:,:,2:end-1) = 1 ./ dlt(2:end) .* ( value(:,3:end,:) - value(:,1:end-2,:));
                res(:,:,end) = (1 ./ dlt(end)) .* (value(:,:,end) - value(:,:,end-1));
                res(:,:,1) = (1 ./ dlt(1)) .* (value(:,:,2) - value(:,:,1));
            end
            
        elseif order == 2
            if dimension == 1
%                 res(2:end-1,:,:) = (1 ./ (dlt)^2) .* ( value(3:end,:,:) + value(1:end-2,:,:) - 2 .* value(2:end-1,:,:));
                res(2:end-1,:,:) = (1 ./ dlt(2:end) .* ( value(3:end,:,:) - value(2:end-1,:,:)) - 1 ./ dlt(1:end-1) .* ( value(2:end-1,:,:) - value(1:end-2,:,:))) ./ dlt(1:end-1);
                res(end,:,:) = (1 ./ dlt(end) ^ 2) .* (value(end,:,:) + value(end-2,:,:) - 2 .* value(end-1,:,:));
                res(1,:,:) = (1 ./ dlt(1)^2) .* (value(3,:,:) + value(1,:,:) - 2 .* value(2,:,:));
            elseif dimension == 2
                res(:,2:end-1,:) = (1 ./ dlt(2:end) .* ( value(:,3:end,:) - value(:,2:end-1,:)) - 1 ./ dlt(1:end-1) .* ( value(:,2:end-1,:) - value(:,1:end-2,:))) ./ dlt(1:end-1);
                res(:,end,:) = (1 ./ dlt(end) ^ 2) .* (value(:,end,:) + value(:,end-2,:) - 2 .* value(:,end-1,:));
                res(:,1,:) = (1 ./ dlt(1)^2) .* (value(:,3,:) + value(:,1,:) - 2 .* value(:,2,:));
            else
                res(:,:,2:end-1) = (1 ./ dlt(2:end) .* ( value(:,:,3:end) - value(:,:,2:end-1)) - 1 ./ dlt(1:end-1) .* ( value(:,:,2:end-1) - value(:,:,1:end-2))) ./ dlt(1:end-1);
                res(:,:,end) = (1 ./ dlt(end) ^ 2) .* (value(:,:,end) + value(:,:,end-2) - 2 .* value(:,:,end-1));
                res(:,:,1) = (1 ./ dlt(1)^2) .* (value(:,:,3) + value(:,:,1) - 2 .* value(:,:,2));
            end
        else
            disp('order not supported');
        end    
        
    end
    if ~isnan(cap)
        res(res < cap) = cap;
    end
    
    
end


