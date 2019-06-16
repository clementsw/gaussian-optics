function C = displace(C,alpha)

%This function displaces a given number of modes. Input the original
%covariance matrix and the complex displacement parameters for all the modes, as a vector.

%Note: the displacement is done in the x and p basis, so if you want to
%displace by alpha, make sure to displace by sqrt(2)*alpha in this code

dim = size(alpha,2);

D = zeros(1,dim*2);
for ii  = 1:dim
   D(ii) =  imag(alpha(ii));
   D(ii+dim) = real(alpha(ii));
end

C.d = C.d + D.';

end