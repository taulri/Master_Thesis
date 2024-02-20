
function [y,bs] = buffered_stats(data,frame,overlap,type)
%   [y] = buffer_variance(data,frame,overlap,type)
%   input  data: input vector of data
%         frame: frame size (scalar value)
%       overlap: overlap size between frames (scalar value)
%          type: i.e. std, mean, median, ...
    if nargin<4
        type = 'std';
    end
    s = buffer(data,frame,overlap);
    bs = size(s,2); % number of segments 
    switch type
        case 'std'
               y = std(s);
        case 'mean'
               y = mean(s);
        case 'median'
               y = median(s);
        case 'var'
               y = var(s);
    end
end

