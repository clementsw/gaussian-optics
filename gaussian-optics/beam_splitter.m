function C = beam_splitter(C,coord,theta)

%This function implements a symmetric beam splitter on two modes, numbered 
%coord(1) and coord(2) with coord(1)<coord(2). The angle of the beam
%splitter is defined by theta.

dim = size(C.M,2)/2;

BS = eye(dim*2);

x = coord(1);
y = coord(2);

BS(x,x) = cos(theta);
BS(x+dim,x+dim) = cos(theta);
BS(y,y) = cos(theta);
BS(y+dim,y+dim) = cos(theta);
BS(x,y+dim) = sin(theta);
BS(y,x+dim) = sin(theta);
BS(x+dim,y) = -sin(theta);
BS(y+dim,x) = -sin(theta);


C.M = BS*C.M*BS.';
C.d = BS*C.d;

end