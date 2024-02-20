function [match] = overlay(ChA, LocA, ChB, LocB, snipsize)
% Get the overlays of events in list A and B, which events do overlay in
% time? Overlap is half the snippet size.

    match = zeros(length(ChB),1);
    indexlist = [];
    for q = 1:length(ChA)

       indexCh = find(ChB == ChA(q));

       spot = (LocB(indexCh) > LocA(q)-(snipsize/2)) + (LocB(indexCh) < LocA(q)+(snipsize/2));
       
       indexChandLoc = find(spot == 2);

       indexlist = [indexlist, indexCh(indexChandLoc)'];
    
    end
    
    match(indexlist) = 1;

end

