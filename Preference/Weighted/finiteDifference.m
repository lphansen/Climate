function res = finiteDifference(value, dimension, order, dlt, scheme, capbound)
% 
% if sum(size(dlt))>2
%     s = size(dlt);
%     if s(1,2)>1
%         dlt=dlt';
%     end
% end
sz = size(value);
res = zeros(sz);
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
    % fetch all step sizes changing points
    points = find(dlt(2:end) - dlt(1:end-1) ~= 0);

    if strcmp(scheme, 'central')
        if order == 1   % first order derivatives
            if dimension == 1
                dlt_multi = reshape(dlt,[length(dlt),1,1]);
                res(2:end-1,:,:) = (1 ./ (dlt_multi(2:end) + dlt_multi(1:end-1))) .* ( value(3:end,:,:) - value(1:end-2,:,:));
                res(end,:,:) = (1 ./ dlt(end)) .* (value(end,:,:) - value(end-1,:,:));
                res(1,:,:) = (1 ./ dlt(1)) .* (value(2,:,:) - value(1,:,:));
                if isempty(points) ~= 0  % modify derivatives at cutoff point
                    for iter = 1 : length(points)
                        pt = points(iter) + 1;
                        left = dlt(pt - 1);
                        right = dlt(pt);
                        if right > left   % left side is the finer grid, build a fake point in the right side
                            vfake = value(pt,:,:) + (value(pt + 1,:,:) - value(pt,:,:)) .* (left / right);
                            res(pt,:,:) = (1 ./ (left * 2)) .* (vfake - value(pt-1,:,:));
                        else % right is the finer grid
                            vfake = value(pt,:,:) - (value(pt,:,:) - value(pt-1,:,:)) .* (right / left);
                            res(pt,:,:) = (1 ./ (right * 2)) .* (value(pt+1,:,:) - vfake);
                        end
                    end
                end
            elseif dimension == 2
                dlt_multi = reshape(dlt,[1,length(dlt),1]);
                res(:,2:end-1,:) = (1 ./ (dlt_multi(2:end) + dlt_multi(1:end-1))) .* ( value(:,3:end,:) - value(:,1:end-2,:));
                res(:,end,:) = (1 ./ dlt(end)) .* (value(:,end,:) - value(:,end-1,:));
                res(:,1,:) = (1 ./ dlt(1)) .* (value(:,2,:) - value(:,1,:));
                if isempty(points) ~= 0  % modify derivatives at cutoff point
                    for iter = 1 : length(points)
                        pt = points(iter) + 1;
                        left = dlt(pt - 1);
                        right = dlt(pt);
                        if right > left   % left side is the finer grid, build a fake point in the right side
                            vfake = value(:,pt,:) + (value(:,pt + 1,:) - value(:,pt,:)) .* (left / right);
                            res(:,pt,:) = (1 ./ (left * 2)) .* (vfake - value(:,pt - 1,:));
                        else % right is the finer grid
                            vfake = value(:,pt,:) - (value(:,pt,:) - value(:,pt - 1,:)) .* (right / left);
                            res(:,pt,:) = (1 ./ (right * 2)) .* (value(:,pt+1,:) - vfake);
                        end
                    end
                end
            else
                dlt_multi = reshape(dlt,[1,1,length(dlt)]);
                res(:,:,2:end-1) = (1 ./ (dlt_multi(2:end) + dlt_multi(1:end-1))) .* ( value(:,:,3:end) - value(:,:,1:end-2));
                res(:,:,end) = (1 ./ dlt(end)) .* (value(:,:,end) - value(:,:,end-1));
                res(:,:,1) = (1 ./ dlt(1)) .* (value(:,:,2) - value(:,:,1));
                if isempty(points) ~= 0  % modify derivatives at cutoff point
                    for iter = 1 : length(points)
                        pt = points(iter) + 1;
                        left = dlt(pt - 1);
                        right = dlt(pt);
                        if right > left   % left side is the finer grid, build a fake point in the right side
                            vfake = value(:,:,pt) + (value(:,:,pt + 1) - value(:,:,pt)) .* (left / right);
                            res(:,:,pt) = (1 ./ (left * 2)) .* (vfake - value(:,:,pt-1));
                        else % right is the finer grid
                            vfake = value(:,:,pt) - (value(:,:,pt) - value(:,:,pt-1)) .* (right / left);
                            res(:,:,pt) = (1 ./ (right * 2)) .* (value(:,:,pt+1) - vfake);
                        end
                    end
                end
            end
            
        elseif order == 2
            if dimension == 1
                dlt_multi = reshape(dlt,[length(dlt),1,1]);
                res(2:end-1,:,:) = ((1 ./ dlt_multi(2:end)) .* ( value(3:end,:,:) - value(2:end-1,:,:)) - (1 ./ dlt_multi(1:end-1)) .* ( value(2:end-1,:,:) - value(1:end-2,:,:))) ./ (0.5 .* (dlt_multi(1:end-1) + dlt_multi(2:end)));
                res(end,:,:) = (1 ./ dlt(end) ^ 2) .* (value(end,:,:) + value(end-2,:,:) - 2 .* value(end-1,:,:));
                res(1,:,:) = (1 ./ dlt(1)^2) .* (value(3,:,:) + value(1,:,:) - 2 .* value(2,:,:));
                if ~isempty(points)  % modify derivatives at cutoff point
                    for iter = 1 : length(points)
                        pt = points(iter) + 1;
                        left = dlt(pt - 1);
                        right = dlt(pt);
                        if right > left   % left side is the finer grid, build a fake point in the right side
                            vfake = value(pt,:,:) + (value(pt + 1,:,:) - value(pt,:,:)) .* (left / right);
                            res(pt,:,:) = (1 ./ (left .^ 2)) .* (vfake + value(pt-1,:,:) - 2 .* value(pt,:,:));
                        else % right is the finer grid, fake point on the left side
                            vfake = value(pt,:,:) - (value(pt,:,:) - value(pt-1,:,:)) .* (right / left);
                            res(pt,:,:) = (1 ./ (right .^ 2)) .* (vfake + value(pt+1,:,:) - 2 .* value(pt,:,:));
                        end
                    end
                end
                
            elseif dimension == 2
                dlt_multi = reshape(dlt,[1,length(dlt),1]);
                res(:,2:end-1,:) = (1 ./ dlt_multi(2:end) .* ( value(:,3:end,:) - value(:,2:end-1,:)) - 1 ./ dlt_multi(1:end-1) .* ( value(:,2:end-1,:) - value(:,1:end-2,:))) ./ (0.5 .* (dlt_multi(1:end-1) + dlt_multi(2:end)));
                res(:,end,:) = (1 ./ dlt(end) ^ 2) .* (value(:,end,:) + value(:,end-2,:) - 2 .* value(:,end-1,:));
                res(:,1,:) = (1 ./ dlt(1)^2) .* (value(:,3,:) + value(:,1,:) - 2 .* value(:,2,:));
                if isempty(points) ~= 0  % modify derivatives at cutoff point
                    for iter = 1 : length(points)
                        pt = points(iter) + 1;
                        left = dlt(pt - 1);
                        right = dlt(pt);
                        if right > left   % left side is the finer grid, build a fake point in the right side
                            vfake = value(:,pt,:) + (value(:,pt + 1,:) - value(:,pt,:)) .* (left / right);
                            res(:,pt,:) = (1 ./ (left .^ 2)) .* (vfake + value(:,pt - 1,:) - 2 .* value(:,pt,:));
                        else % right is the finer grid
                            vfake = value(:,pt,:) - (value(:,pt,:) - value(:,pt - 1,:)) .* (right / left);
                            res(:,pt,:) = (1 ./ (right .^ 2)) .* (value(:,pt+1,:) + vfake - 2 .* value(:,pt,:));
                        end
                    end
                end
            else
%                 for ite = 2 : length(res(1,1,:)) - 1
%                     res(:,:,ite) = ((value(:,:,ite + 1) - value(:,:,ite)) ./ dlt(ite) - (value(:,:,ite) - value(:,:,ite-1))./ dlt(ite-1))./ ( 0.5*(dlt(ite) + dlt(ite-1)));
%                 end
                dlt_multi = reshape(dlt,[1,1,length(dlt)]);
                res(:,:,2:end-1) = (1 ./ dlt_multi(2:end) .* ( value(:,:,3:end) - value(:,:,2:end-1)) - 1 ./ dlt_multi(1:end-1) .* ( value(:,:,2:end-1) - value(:,:,1:end-2))) ./ (0.5 .* (dlt_multi(1:end-1) + dlt_multi(2:end)));
                res(:,:,end) = (1 ./ dlt(end) ^ 2) .* (value(:,:,end) + value(:,:,end-2) - 2 .* value(:,:,end-1));
                res(:,:,1) = (1 ./ dlt(1)^2) .* (value(:,:,3) + value(:,:,1) - 2 .* value(:,:,2));
                if isempty(points) ~= 0  % modify derivatives at cutoff point
                    for iter = 1 : length(points)
                        pt = points(iter) + 1;
                        left = dlt(pt - 1);
                        right = dlt(pt);
                        if right > left   % left side is the finer grid, build a fake point in the right side
                            vfake = value(:,:,pt) + (value(:,:,pt + 1) - value(:,:,pt)) .* (left / right);
                            res(:,:,pt) = (1 ./ (left .^ 2)) .* (vfake + value(:,:,pt-1) - 2 .* value(:,:,pt));
                        else % right is the finer grid
                            vfake = value(:,:,pt) - (value(:,:,pt) - value(:,:,pt-1)) .* (right / left);
                            res(:,:,pt) = (1 ./ (right .^ 2)) .* (value(:,:,pt+1) + vfake - 2 .*value(:,:,pt));
                        end
                    end
                end
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


