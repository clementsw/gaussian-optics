%Yields the photon number statistics of a Gaussian state with covariance
%matrix C.M and first order moments C.d. resolution is the highest number of
%photons to be resolved. 

%The following can only really be understood if read in conjunction with
%"Multidimensional Hermite polynomials and photon distributions for
%polymode mixed light", by Dodonov et al, PRA 1994

%Note: this computation 1) scales exponentially with the number of modes
%and 2) is not performed in this code in a particularly efficient way.
%Calculating meaningful photon statistics for more than about 4 modes can
%take a long time

function P = photon_number_stats(C,resolution)

M = C.M;
d = C.d;

size_d = size(M);
dim = size_d(1); %dimension of the problem

%In the paper
U=eye(dim);
for ii = 1:dim/2
    U(ii,ii) = -1i;
    U(ii+dim/2,ii) = 1;
    U(ii,ii+dim/2) = 1i;
end
U=U/sqrt(2);

U1 = U';
U2 = conj(U);
U3 = U.';

%In the paper
R=U1*(eye(dim)-2*M)*inv(eye(dim)+2*M)*U2;
y = 2*U3*inv(eye(dim)-2*M)*d;

%Create dim-dimensional tensor in which every entry is the value of the corresponding
%Hermite polynomial at position y
X= zeros(1,dim);
for ii = 1:dim
    X(ii) = resolution;
end
H=zeros(X);

%first element, H(0,0,0,0)
x=num2cell(ones(1,dim));
H(x{:}) = 1;
 
%In the following, I go through the H tensor and fill it in element by
%element, by "jumping" from element n to element n+ek (see recursion
%relation to understand how this works)

%Initialise
nextPos = ones(1,dim);
jumpFrom = ones(1,dim);
jump = 0;

for jj = 1:resolution^dim - 1
    
    %Figure out what the next position is and which position to jump from
    %to get there
    if jump > 0
        jumpFrom(jump+1) = jumpFrom(jump+1) + 1;
        jump = 0;
    end
    
    for ii = 1:dim
        forwardStep=zeros(1,dim);
        forwardStep(ii) = 1;
        if forwardStep(ii) + nextPos(ii) > resolution
            nextPos(ii) = 1;
            jumpFrom(ii) = 1;
            jump = ii;
        else
            jumpFrom(ii) = nextPos(ii);
            nextPos(ii) = nextPos(ii) + 1;       
            break
        end
    end
    
    
    %Implement recursion relation
    ek = nextPos-jumpFrom;
    k = find(ek);
    nextCoordinate = num2cell(nextPos);
    fromCoordinate = num2cell(jumpFrom);
    
    for ii = 1:dim 
        H(nextCoordinate{:}) = H(nextCoordinate{:}) + R(k,ii)*y(ii);
    end
    H(nextCoordinate{:}) = H(nextCoordinate{:})*H(fromCoordinate{:});
    
    for ii = 1:dim 
        
        if jumpFrom(ii) > 1
            prevJump = zeros(1,dim);
            prevJump(ii) = 1;
            prevCoordinate = num2cell(jumpFrom - prevJump);
            H(nextCoordinate{:}) = H(nextCoordinate{:}) - (jumpFrom(ii)-1)*R(k,ii)*H(prevCoordinate{:});
        end
    end
end

H=real(H);

%Now calculate the probabilities

%Calculate a tensor of probabilities
X= zeros(1,dim/2);
for ii = 1:dim/2
    X(ii) = resolution;
end
P=zeros(X);

%Annoying syntax detail if dim/2==1
if dim == 2
    P = zeros(1,resolution);
end

%first element, P(0,0,0,...). See paper.
x=num2cell(ones(1,dim/2));
P0 = exp(-d'*inv(2*M+eye(dim))*d)/sqrt(det(M+0.5*eye(dim)));
P(x{:}) = P0;

%Similar procedure to what precedes. Move forward in the P tensor and fill
%it element by element.
nextP = ones(1,dim/2);
for jj = 1:resolution^(dim/2) - 1
    
    %Figure out what the next coordinate to fill in is
    for ii = 1:dim/2
        jumpTo=zeros(1,dim/2);
        jumpTo(ii) = 1;
        if nextP(ii) + jumpTo(ii) > resolution
            nextP(ii) = 1;
        else
            nextP(ii) = nextP(ii) + 1;
            break
        end
    end
    
    %Calculate probability. See paper for the formula.
    nextCoord = num2cell(nextP);
    whichH = zeros(1,dim);
    for kk = 1:dim/2
        whichH(kk) = nextP(kk);
        whichH(kk+dim/2) = nextP(kk);
    end
    whichH = num2cell(whichH);
    P(nextCoord{:}) = P0*H(whichH{:});
    for  kk = 1:dim/2
        P(nextCoord{:}) = P(nextCoord{:})/factorial(nextP(kk)-1);
    end
    
end

end