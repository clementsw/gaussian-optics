function C = add_loss(C,t)

%This function adds loss to a given number of modes. Input the original
%covariance matrix and the losses (t=0: no loss; t=1: no transmission) for all the modes, as a vector.

dim = size(t,2);

loss = ones(1,2*dim);

for ii = 1:dim
    loss(ii) = sqrt(1-t(ii));
    loss(ii+dim) = sqrt(1-t(ii));
end

T=diag(loss);

C.M = T*C.M*T + (eye(2*dim)-T^2)/2;
C.d = T*C.d;

end