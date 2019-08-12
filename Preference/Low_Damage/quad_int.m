function result = quad_int(f, a, b, n, method)
%This function takes a function f to integrate from the multidimensional
%interval specified by the row vectors a and b. N different points are used
%in the quadrature method. Legendre and Hermite basis functions are
%currently supported. In the case of Hermite methodology b is the normal
%density and a is the normal mean.

% Created by John Wilson (johnrwilson@uchicago.edu)

if strcmp(method,'legendre')==1
    [xs, ws] = quad_points_legendre(n);
    g = @(x) f((b - a) / 2 .* x + (a+b)/2);
    s = prod((b - a) ./ 2);
    
elseif strcmp(method,'hermite')==1
    [xs, ws] = quad_points_hermite(n);
    g = @(x) f(sqrt(2) * b * x + a);
    s = 1 / sqrt(pi);
    
else
    ME = MException('ValueError:BadParameter', ...
        'Invalid method parameter used in quadrature.');
    throw(ME);
end

dim = size(a, 2);
sum = 0;
if dim == 3
    for i = 1:n
        for j = 1:n
            for k = 1:n
                sum = sum + (ws(i) * ws(j) * ws(k)) .* g([xs(i) xs(j) xs(k)]);
            end
        end
    end
end

if dim == 2
    for i = 1:n
        for j = 1:n
            sum = sum + (ws(i) * ws(j)) .* g([xs(i) xs(j)]);
        end
    end
end

if dim == 1
    for i = 1:n
        sum = sum + ws(i) .* g(xs(i));
    end
end

result = s * sum;
end