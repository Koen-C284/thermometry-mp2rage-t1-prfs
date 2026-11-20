function A = get_zernike_function(dims, order, safe)
    if nargin<3, safe = false; end
    [yy,xx] = meshgrid(linspace(-1,1,dims(1)), linspace(-1,1,dims(2)));
    x = xx(:); y = yy(:);
    [phi, r] = cart2pol(x,y);
    r = min(r,1);

    cols = {ones(numel(r),1)};  % Z_0^0
    for n = 1:order
        for m = -n:2:n          % Zernike parity: n-m even
            if (mod(n-m,2)~=0), continue; end
            if safe && m==0     % drop axis-symmetric terms -> no center bowls
                continue
            end
            Rnm = radialZernike(n,abs(m),r);
            if m>0
                cols{end+1} = Rnm.*cos(m*phi); %#ok<AGROW>
            elseif m<0
                cols{end+1} = Rnm.*sin(-m*phi); %#ok<AGROW>
            end
        end
    end
    A = cat(2, cols{:});
end

function R = radialZernike(n,m,r)
    R = zeros(size(r));
    for s = 0:((n-m)/2)
        c = (-1)^s * factorial(n-s) / ...
            ( factorial(s)*factorial((n+m)/2 - s)*factorial((n-m)/2 - s) );
        R = R + c * r.^(n-2*s);
    end
end
