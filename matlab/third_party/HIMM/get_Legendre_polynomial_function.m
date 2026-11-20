function A = get_Legendre_polynomial_function(dims, order, safe)
% Safe 2D Legendre tensor basis up to "order".
% Excludes P_{2k}(x)P_0(y) and P_0(x)P_{2k}(y) when safe=true.
    if nargin < 3, safe = false; end
    [yc,xc] = meshgrid(linspace(-1,1,dims(1)), linspace(-1,1,dims(2)));
    xc = xc(:); yc = yc(:);

    Px = legendre_polys_1d(xc, order);
    Py = legendre_polys_1d(yc, order);

    cols = {};
    for nx = 0:order
        for ny = 0:(order-nx)
            if safe && ( (nx>0 && ny==0 && mod(nx,2)==0) || (ny>0 && nx==0 && mod(ny,2)==0) )
                continue; % exclude pure even in one axis
            end
            cols{end+1} = Px(:,nx+1).*Py(:,ny+1); %#ok<AGROW>
        end
    end
    A = cat(2, cols{:});
end

function P = legendre_polys_1d(x, order)
% Evaluate Legendre P_0..P_order at points x (column vector).
    N = numel(x);
    P = zeros(N, order+1);
    P(:,1) = 1;
    if order >= 1, P(:,2) = x; end
    for n = 1:order-1
        P(:,n+2) = ((2*n+1).*x.*P(:,n+1) - n*P(:,n))/(n+1);
    end
end