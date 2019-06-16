function C = amplify(C,x)

%This function amplifies a given number of modes in a phase insensitive
%way, by doing two mode squeezing with squeezing parameter r on the target
%mode and on an ancilla vacuum, and by tracing over the vacuum. Input the original
%covariance matrix and the "squeezing parameters" for all the modes.

dim = size(x,2);

D = zeros(1,dim*2);
for ii  = 1:dim
   D(ii) =  cosh(x(ii));
   D(ii+dim) = cosh(x(ii));
end
D = diag(D);

S = zeros(1,dim*2);
for ii  = 1:dim
   S(ii) =  sinh(x(ii));
   S(ii+dim) = sinh(x(ii));
end
S = diag(S);

C.M = D*C.M*D.' + S^2/2;
C.d = D*C.d;

end