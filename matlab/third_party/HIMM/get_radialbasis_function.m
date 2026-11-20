function A = get_radialbasis_function(dims, order, safe)
% Gaussian RBF grid; "order" controls grid density (order+1 per axis).
% Excludes center RBF when safe=true.
    if nargin < 3, safe = false; end
    [yc,xc] = meshgrid(linspace(-1,1,dims(1)), linspace(-1,1,dims(2)));
    xc = xc(:); yc = yc(:);

    nper = max(1, order+1);
    gmin = -0.7; gmax = 0.7;
    cgrid = linspace(gmin,gmax,nper);
    [cy,cx] = meshgrid(cgrid,cgrid);
    centers = [cx(:), cy(:)];

    if safe
        centers(all(abs(centers) < 1e-12, 2),:) = []; % drop (0,0)
    end

    if numel(cgrid) > 1
        sigma = 0.6 * mean(diff(cgrid));
    else
        sigma = 0.5; % fallback
    end

    cols = {ones(numel(xc),1)}; % include constant
    for k = 1:size(centers,1)
        r2 = (xc-centers(k,1)).^2 + (yc-centers(k,2)).^2;
        cols{end+1} = exp(-r2/(2*sigma^2)); %#ok<AGROW>
    end
    A = cat(2, cols{:});
end