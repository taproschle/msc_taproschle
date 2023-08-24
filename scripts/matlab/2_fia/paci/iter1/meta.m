function dxdt   = meta(ModelName,t,x,theta,dim)

%This is an auxiliary function for simulation using the provided Zenteno model in symbolical format.

%Dont change these
nx  =   dim(1); 
nth =   dim(2);

%Model states
X   = x(1:nx);

%Construction of d(x_theta)/dt from eq. 3 (Stigter & Molenaar, 2015) (Dont change unless inputs (Ui) must be modified)
f=eval(['@' ModelName]);
xdot=f(t,X,theta);

dxdth=reshape(x((nx+1):end),nx,nth);

dfdxModel=eval(['@dfdx_' ModelName]);
dfdx=dfdxModel(t,X,theta);

dfdthModel=eval(['@dfdth_' ModelName]);
dfdth=dfdthModel(t,X,theta);

dxdthdot=dfdx*dxdth+dfdth;

dxdt=[xdot; dxdthdot(:)];