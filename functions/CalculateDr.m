function Dr = CalculateDr(A,delta,long,npointsMSDset)
% INPUT
% A = x y framenumber
% delta = dt
% long = length track
% npointsMSDset = smaller than minimal track length, sliding window for MSD is (npointsMSDset + 1) frames, default 10

S = size(A,1);

% calculate instantaneous diffusion coefficient. Middle point in window
npointsMSD = min(long-1,npointsMSDset);
hp = round(0.5 * npointsMSD);
for j=1:S
    if j<=hp
        Dins(j)=MSqDin(A(1:1+npointsMSD,:),delta,4);
    elseif j>S-hp 
        Dins(j)=MSqDin(A(S-npointsMSD:S,:),delta,4);
    else
        Dins(j)=MSqDin(A(j-hp:j-hp+npointsMSD,:),delta,4);
    end
    if Dins(j)<0
        Dins(j)=0;
    end
end

Dr = Dins;