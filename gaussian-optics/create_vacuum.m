function C = create_vacuum(N)

%Create the covariance matrix and displacement vector for an N-mode vacuum.

C.M = eye(2*N)/2;
C.d = zeros(2*N,1);

%photonNumberStats goes crazy if there is no squeezing at all. The
%following solves that problem
S = ones(1,N)*0.000001;
C=squeeze(C,S);

end