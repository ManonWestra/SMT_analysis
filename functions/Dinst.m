function D = Dinst(time, msd, np)
    P = polyfit(time(1:np), msd(1:np), 1);   
    if P(1) <0
        D = [];
    else
        D = P(1) / 4;
    end
    