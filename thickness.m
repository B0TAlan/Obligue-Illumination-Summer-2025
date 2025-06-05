function THICC = thickness(normComplexFeild, lambda, ns, nim, check)
arguments
    normComplexFeild 
    lambda 
    ns 
    nim 
    check
end
phaseDeg = angle(normComplexFeild);
[N,M] = size(normComplexFeild);

THICC= (phaseDeg * lambda)/ (2* pi * (ns-nim));

if check == 1 % creates plot to check if thickness center allines w/ center
    x = [M/2,M];% creates the range for the ref line
    y = [N/2, N/2];

    prefil = improfile(THICC, x,y);
    
    figure(8), imagesc(THICC), colormap hot, title '0 deg thickness'
    hold on
    plot(x, y, 'c-', 'LineWidth', 2); % 'c-' for cyan line
    hold off
    

    %plot of the profiles
    figure(6)
    plot(prefil)
elseif check == 2
    x = [M/2,M];% creates the range for the ref line
    y = [N/2+1, N/2];

    prefil = improfile(THICC, x,y);
    
    figure(8), imagesc(THICC), colormap hot, title '0 deg thickness'
    hold on
    plot(x, y, 'c-', 'LineWidth', 2); % 'c-' for cyan line
    hold off
    

    %plot of the profiles
    figure(6)
    plot(prefil)
end
end