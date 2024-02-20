function [inclusion] = select_candidates_HFO(error, V, per, pos, maxV)
% Takes the approximation error and v-factor as input and decides whether snippet is
% candidate for localization of epileptogenic zone

inclusion = 0;
error = 1-error; 
slopes = diff(error);

if all(slopes >= 0) % Is the slope always above = 0? Should be!    
        if error(pos) > per % Setting a threshold for the reconstruction error
            if error(length(error)) > 0.8 % Setting threshold for the final reconstruction error 
                if max(V) < maxV % add constrains on V-factor
                    inclusion = 1;
                end
            end
        end 
end
