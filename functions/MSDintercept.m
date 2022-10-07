function [MSDslope, D,  MSDint] = MSDintercept(time, msd, np)
    P = polyfit(time(2:np), msd(2:np), 1);   
    MSDslope = P(1);
    MSDint = P(2);
    if P(1) <0
        D = 100;
    else
        D = P(1) / 4;
    end
   % Dlog = log10(D);
    end
    