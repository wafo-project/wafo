function [x,y] = wdiscretize(fun, a, b, varargin)
%    Automatic discretization of function
%
%    Parameters
%    ----------
%    fun : callable
%        function to discretize
%    a,b : real scalars
%        evaluation limits
%    tol : real, scalar, default 0.005
%        absoute error tolerance
%    n : scalar integer, default n=5
%        number of values to start the discretization with.
%    method : string, default 'linear'
%        defining method of gridding, options are 'linear' and 'adaptive'
%    nmax : scalar integer default 2 ** 15 + 1
%        maximum number of values 
%    maxtries : scalar integer default 5
%        maximum number of attempts to decrease the maximum error before
%        termination.
%
%    Returns
%    -------
%    x : discretized values
%    y : fun(x)
%
% Example
%  [x,y] = wdiscretize(@(x)cos(x), 0, pi);
%  [xa, ya] = wdiscretize(@(x)cos(x), 0, pi, 'method', 'adaptive');
%  assert(xa(1:5), ... 
%      [ 0.        ,  0.19634954,  0.39269908,  0.58904862,  0.78539816], 1e-7);
%
%  plot(x, y, xa, ya, 'r.');
% 
%  close('all');
    options = struct('tol', 0.005, 'n', 5, 'method', 'linear', 'nmax', 2^15+1, 'maxtries', 5);
 
    opt = parseoptions(options, varargin{:});
  
    if opt.method(1)=='a',
        [x,y] = discretize_adaptive(fun, a, b, opt.tol, opt.n, opt.nmax, opt.maxtries);
    else
        [x,y] = discretize_linear(fun, a, b, opt.tol, opt.n, opt.nmax, opt.maxtries);
    end
end


function [x, y] = discretize_linear(fun, a, b, tol, n, nmax, maxtries),
    % Automatic discretization of function, linear gridding
    
    x = linspace(a, b, n);
    y = fun(x);

    err0 = inf;
    err = 10000;
    TINY = realmin;
    num_tries = 0;
    while (num_tries<maxtries && err > tol && n < nmax),
        err0 = err;
        x0 = x;
        y0 = y;
        n = 2 * (n - 1) + 1;
        x = linspace(a, b, n);
        y = fun(x);
        y00 = interp1(x0, y0, x, 'linear');
        err = 0.5 * max(abs(y00 - y) ./ (abs(y00) + abs(y) + TINY + tol));
        num_tries = num_tries + (abs (err - err0) <= tol/2);
    end
    return 
end


function [x,fx] = discretize_adaptive(fun, a, b, tol, n, nmax, maxtries),
    % Automatic discretization of function, adaptive gridding.
    
    n = n + (mod(n, 2) == 0);  % make sure n is odd
    x = linspace(a, b, n);
    fx = fun(x);

    n2 = (n - 1) / 2;
    erri_ = [zeros(1, n2); ones(1,n2)];
    erri = erri_(:).';
    err = max(erri);
    err0 = inf;
    TINY = realmin;
    num_tries = 0;
    % while (err != err0 and err > tol and n < nmax):
    for j = [1:50],
        if num_tries<maxtries && err > tol && length(fx) < nmax,
            err0 = err;
            % find top errors

            I = find(erri > tol);
            % double the sample rate in intervals with the most error
            y_ = [(x(I) + x(I - 1)) / 2; (x(I + 1) + x(I)) / 2];
            y = y_(:).';
            fy = fun(y);
            fy0 = interp1(x, fx, y, 'linear');
            erri = 0.5 * (abs(fy0 - fy) ./ (abs(fy0) + abs(fy) + TINY + tol));
            err = max(erri);

            x = [x, y];
            [x,I] = sort(x);
            erri_ = [zeros(1,length(fx)), erri];
            erri = erri_(I);
            fx_ = [fx, fy];
            fx = fx_(I);

            num_tries = num_tries + (abs(err - err0) <= tol/2);
        else,
            break
        end
    end
    if any(erri > tol),
        warning('WAFO:WDISCRETIZE', ...
          sprintf('Adaptive discretization did not converge! Check function for sharp edges! (max error = %g, n=%d) ', err, length(fx)))
    end
    if j>=50,
        warning('WAFO:WDISCRETIZE',...
          sprintf('Recursion level limit reached j=%d', j));
    end
    return
end