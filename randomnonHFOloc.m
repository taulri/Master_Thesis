function [randomloc] = randomnonHFOloc(n, listlocAndCh, size_atom, sizedata, channeln)

aloc = (size_atom/2)+1;
bloc = sizedata -(size_atom/2); 
ach = 1;
bch = channeln;
randomloc = zeros(n,2);

for k = 1:n
    xchR = randi([ach,bch]);
    ind = find(listlocAndCh(:,2) == xchR);
    cand = listlocAndCh(ind,:);
    if size(cand,1)==0
        xlocR = randi([aloc,bloc]);
        randomloc(k,1) = xlocR;
        randomloc(k,2) = xchR;
    elseif size(cand,1)==1
        if cand(1,1)+2000 < sizedata
            xlocR = cand(1,1)+1400;
            randomloc(k,1) = xlocR;
            randomloc(k,2) = xchR;
        else 
            xlocR = cand(1,1)-2000;
            randomloc(k,1) = xlocR;
            randomloc(k,2) = xchR;
        end
    else
        dx = diff(cand(:,1));
        if ~any(dx > 5000)
            if cand(end,1)+2000 < sizedata
                xlocR = cand(end,1)+1400;
                randomloc(k,1) = xlocR;
                randomloc(k,2) = xchR;
            end
        end
        ok = 0;
        while ok == 0
            lendx = length(dx);
            x = randi(lendx);
            if abs(dx(x)) > 5000
                xlocR = cand(x,1)+2500;
                ok = 1;
                randomloc(k,1) = xlocR;
                randomloc(k,2) = xchR;
            end
        end
    end
end

end