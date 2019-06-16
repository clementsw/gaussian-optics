function C = squeeze(C,x)

%This function squeezes a given number of modes. Input the original
%covariance matrix and the squeezing parameters for all the modes.

dim = size(x,2);

D = zeros(1,dim*2);
for ii  = 1:dim
   D(ii) =  exp(-x(ii));
   D(ii+dim) = exp(x(ii));
end

S=diag(D);

C.M = S*C.M*S.';
C.d = S*C.d;

end