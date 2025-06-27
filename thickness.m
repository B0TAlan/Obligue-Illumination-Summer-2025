
% rewrite so that input is phaseDeg 
% make it more universal 


function THICC = thickness(phaseDeg,x,y, lambda, ns, nim, check)
arguments
    phaseDeg
    x {mustBeNumeric}
    y {mustBeNumeric}
    lambda {mustBeNumeric}
    ns {mustBeNumeric}
    nim {mustBeNumeric}
    check {mustBeInRange(check, 0, 1)}
end

THICC= (phaseDeg * lambda)/ (2* pi * (ns-nim));

if check == 1 % creates plot to check if thickness center allines w/ center
    %disp(1)
    
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