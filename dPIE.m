function [d,d_down,d_up]=dPIE(Obj_T,d,d_down,d_up,L0,lamda,propagator,delta_L0)
    delta_z=lamda*(2*d/L0)^2;
    K=10;
    p=0.01;
    U1=0;U2=0;
    i=1;
    for k=-K/2:K/2
        if k~=0
            d1=k*delta_z;
            Ud_sm=Propagate(Obj_T,propagator,delta_L0,lamda,-d1);
            [Nx,Ny]=size(Ud_sm);
            dx=zeros(Nx,Ny);
            dy=zeros(Nx,Ny);
            dy(1:Nx-1,1:Ny)=(Ud_sm(2:Nx,1:Ny)-Ud_sm(1:Nx-1,1:Ny));
            dx(1:Nx,1:Ny-1)=(Ud_sm(1:Nx,2:Ny)-Ud_sm(1:Nx,1:Ny-1));
            U=sqrt((abs(dx)).^2+(abs(dy)).^2+p);
            u(i)=sum(sum(U(:)));
            U1=U1+u(i)*k;
            i=i+1;
        end
    end
    sita_z=U1;
    if sita_z>0
        d_up=d;
    else
        d_down=d;
    end
    d=(d_down+d_up)/2;
end