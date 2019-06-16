% Calculate Wigner functions for single mode gaussian states, as a function
% of the covariance matrix and vector of first order moments

function W = calculate_single_mode_wigner(C,size,step)

    M = C.M;
    d = C.d;
    
    %Wigner function
    W=zeros(size,size);
    
    for ii = 1:size
        for jj=1:size
            x=step*(jj-size/2);
            q=step*(ii-size/2);
            
            z = [x-d(1),q-d(2)]*inv(M)*[x-d(1),q-d(2)].';
            W(ii,jj) = exp(-z/2)/(2*pi*sqrt(det(M)));
            
        end
    end

end
