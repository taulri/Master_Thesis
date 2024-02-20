function [cross] = crossing(res, min_freq, fs, Nc,th)

N = round(1/min_freq*fs); % Minimal distance so that it's not the same event 
h = ones(Nc-1,1); % Array with minimal number of crossings as ones 
index_HFO = [];
[~,~,ix] = zerocrossrate(res-th'); % logical of positions where we have a crossing
index_cross = find(ix); % find indicies of crossing 
if ~isempty(index_cross) % if there are crossings
    dist_cross = diff(index_cross); % distance of indicies of crossing 
    if ~isempty(dist_cross)
        x_HFO = conv(dist_cross<N,h); % makes sure distance is smaller than 25 within HFO
        index_HFO = find(x_HFO(5:end) == Nc-1); % index of HFO in crossing array 
        index_HFO_final = index_cross(index_HFO); % time index of HFO in filtered signal       
        % Check that the indexes are not the same event
        if length(index_HFO_final)>= 2
            index_HFO_final = index_HFO_final([1,find(diff(index_HFO_final) > N)+1]);
        end
        n_HFO = length(index_HFO_final); % number of detected HFO
    else
        n_HFO = 0;
    end
else 
    n_HFO = 0;
end
cross = n_HFO<1;

