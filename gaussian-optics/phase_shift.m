function C = phase_shift(C,phi)

%This function adds a phase shift to a given number of modes. Input the original
%covariance matrix and the phase shifts for all the modes, as a vector.

dim = size(phi,2);

R=eye(2*dim);

for ii = 1:dim
    R(ii,ii) = cos(phi(ii));
    R(ii,dim+ii) = -sin(phi(ii));
    R(dim+ii,ii) = sin(phi(ii));
    R(dim+ii,dim+ii) = cos(phi(ii));
end
    
C.M = R*C.M*R.';
C.d = R*C.d;

end