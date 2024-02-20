function [y,coeff,loc,residual,Error]=OMP_process(D,X,L,S1, fs, minHz, Nc, th) 
% Sparse coding of the signal based on a given 
% dictionary and specified number of atoms to use. 
% input arguments: 
%               D - the dictionary (its columns MUST be normalized)
%               X - the signal to represent
%               L - the max. number of iterations for each signal -> 10
%               S1 - the min Error to stop OMP -> 0.05
%               fs - sampling frequency
%               minHz - lower bound of HFO range 
%               th - initial threshold for snippet 


[~,P]=size(X);
[~,K]=size(D);
Error_Diff = [];
for k=1:1:P % We have only one signal to represent 
    a=[];
    x=X(:,k); % Chose this singal 
    residual(:,1)=x; % First residual is this signal 
    indx = [];
    for j=1:1:L % iteration max is 10
        proj=D'*residual(:,j); % calculate correlation of signal with a dictionary atom 
        [~,pos]=max(abs(proj)); % take the max correlation 
        pos=pos(1);
        indx=[indx pos];
        a=pinv(D(:,indx(1:j)))*x; % calculate the coefficient 
        residual(:,j+1)=x-D(:,indx(1:j))*a;
        y(:,j) = D(:,indx(1:j))*a;
        %Error       
        Error(j)=sqrt(sum(residual(:,j+1).^2))/sqrt(sum(X.^2)); 
        temp=zeros(K,1);
        temp(indx(1:j))=a;
        coeff(:,k)=(temp);
        loc = indx;
        % Error Breakpoint
        if (Error(j) < S1)
            break;
        end
        if crossing(residual(:,j+1), minHz, fs, Nc, th) % Check if the residual still crosses the threshold 
            break;
        end
       
    end 
end
return;