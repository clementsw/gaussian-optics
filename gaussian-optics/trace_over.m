function C=trace_over(C,ii) 

%Traces over all the modes, keeping only mode ii. The Wigner function of
%the remaining mode can then be used with calculate_single_mode_wigner

s=size(C.M,2);

M = zeros(2,2);
d = zeros(2,1);

M(1,1) = C.M(ii,ii);
M(1,2) = C.M(ii,s/2+ii);
M(2,1) = C.M(s/2+ii,ii);
M(2,2) = C.M(s/2+ii,s/2+ii);

d(1) = C.d(ii);
d(2) = C.d(s/2+ii);

C.M = M;
C.d = d;
end