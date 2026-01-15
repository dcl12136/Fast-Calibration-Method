function [Obj_T,illumR0,Obj_Fun,Error,RMS_AM,RMS_PH]=Ptychography(propagator,iteration,len,Obj_S0,Obj_T,illumR0,Obj_Fun,lamda,delta_L0,Dect_Std,illum_lim,d,x1,y1,method,sita,alpha,beta,Object0)
% propagator: one of 'fourier', 'fresnel' or 'angular spectrum'
% N: pixel number
% iteration: number of the seting iteration
% len : number of scaning positions
% Obj_S0=Rp*2+(Num-1)*Obj_Mov; the range of the tested object
% Obj_T and illumR0: intial guesess of the object and probe
% delta_L0: pixel size
% Dect_Std and illum_lim: the recorded intensity with object and without object
% d: the axial distance between the object and CCD
% x1 and y1: the coordinates of the scaning position
% method, (alpha, beta): updating function and updating coefficient
% sita: stopping criterion for itertaion
%parameters for mPIE
% Eta_Obj=0.5;Eta_Pro=0.75;
% V_Obj=0;V_Pro=0;
% Robj=0.2;Rpro=1;
% T=10;
for i=1:iteration
    sum_obj=0; Weight=0;
    for m1=1:len
        %         Obj_Fun=zeros(N,N);
        Obj_Fun(x1(m1):x1(m1)+Obj_S0-1,y1(m1):y1(m1)+Obj_S0-1)=Obj_T;
        Input_L=illumR0.*Obj_Fun;
        Dect_R=Propagate(Input_L,propagator,delta_L0,lamda,d);
        Weight=Weight+sum(sum((abs(Dect_R).^2-Dect_Std(:,:,m1).^2).^2));
        sum_obj=sum_obj+sum(sum(Dect_Std(:,:,m1).^4));
        Dect_R1=Dect_Std(:,:,m1) .*exp(1i*angle(Dect_R));
        Return_L=Propagate(Dect_R1,propagator,delta_L0,lamda,-d);
        switch method     % updating function
            case 'rPIE'
                denomO=(1-alpha)*abs(illumR0.*conj(illumR0))+alpha*max(max(illumR0.*conj(illumR0)));
                Obj_Fun=Obj_Fun+conj(illumR0).*(Return_L-Input_L)./denomO;
                denomP=(1-beta).*abs(Obj_Fun.*conj(Obj_Fun))+beta*max(max(Obj_Fun.*conj(Obj_Fun)));
                illumR0=illumR0+conj(Obj_Fun).*(Return_L-Input_L)./denomP;
                Obj_T=Obj_Fun(x1(m1):x1(m1)+Obj_S0-1,y1(m1):y1(m1)+Obj_S0-1);
            case 'ePIE'
                denomO=max(max(illumR0.*conj(illumR0)))+eps;
                Obj_Fun=Obj_Fun+alpha*conj(illumR0).*(Return_L-Input_L)./denomO;
                denomP=max(max(Obj_Fun.*conj(Obj_Fun)))+eps;
                illumR0=illumR0+beta*conj(Obj_Fun).*(Return_L-Input_L)./denomP;
                Obj_T=Obj_Fun(x1(m1):x1(m1)+Obj_S0-1,y1(m1):y1(m1)+Obj_S0-1);
            case 'mPIE'
                denomO=(1-alpha)*abs(illumR0).^2+alpha*max(abs(illumR0(:))).^2;
                Obj_Fun=Obj_Fun+Robj*conj(illumR0).*(Return_L-Input_L)./denomO;
                denomP=(1-beta).*abs( Obj_Fun).^2+beta*max(abs(Obj_Fun(:))).^2;
                illumR0=illumR0+Rpro*conj( Obj_Fun).*(Return_L-Input_L)./denomP;
                Obj_T=Obj_Fun(x1(m1):x1(m1)+Obj_S0-1,y1(m1):y1(m1)+Obj_S0-1);
                obj_T(:,:,m1)=Obj_T;
                illumR01(:,:,m1)=illumR0;
                if mod(m1,T)==0
                    V_Obj=Eta_Obj*V_Obj+(obj_T(:,:,m1)-obj_T(:,:,m1+1-T));
                    Obj_T=obj_T(:,:,m1)+Eta_Obj*V_Obj;
                    V_Pro=Eta_Pro*V_Pro+(illumR01(:,:,m1)-illumR01(:,:,m1+1-T));
                    illumR0=illumR01(:,:,m1)+Eta_Pro*V_Pro;
                end
                %%%%%%%%%%%%% amplitude restricttion for illumination %%%%%%%%%%%%%
                illum_Tr = Propagate (illumR0,propagator,delta_L0,lamda,d);
                illum_Tr0=illum_lim.*exp(1i*angle(illum_Tr));
                illumT = Propagate (illum_Tr0,propagator,delta_L0,lamda,-d);
                illumR0=illumT;
        end
        Error(i)=(Weight/sum_obj);
        dAM=(abs(Obj_T)-abs(Object0))^2;
        dPH=(angle(Obj_T)-angle(Object0))^2;
        AM=dAM(200:500,200:500);PH=dPH(200:500,200:500);
        [m,n]=size(AM);
        RMS_AM(i)=sqrt(sum(sum(AM))/(m*n));
        RMS_PH(i)=sqrt(sum(sum(PH))/(m*n));
%         if  (i>1 && abs(Error(i)-Error(i-1))/Error(i-1)<sita) || (i>iteration)
%             break;
%         end
    end
end