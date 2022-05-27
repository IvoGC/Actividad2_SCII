clc; clear ; close all
%{
---------------------------------------------------------------------------
Calcular un controlador que haga evolucionar al péndulo en el equilibrio estable, partiendo 
de una condición inicial nula en el desplazamiento y el ángulo en pi que termine en 10 metros evitando 
las oscilaciones de la masa m, considerando que es una grúa. Una vez que delta=10 modificar a m a un 
valor 10 veces mayor y volver al origen evitando oscilaciones.
---------------------------------------------------------------------------
%}
%DEFINO PARAMETROS
m = 0.1; F = 0.1; l = 0.6; g = 9.8; M = 0.5; 

%DEFINO MATRICES
%X=[delta delta_p phi phi_p]
%Sistema linealizado en el punto Xo
X_o=[0 0 0 0];%equilibio inestable phi=0
A = [0,1,0,0 ; 0,(-F/M),(-m*g/M),0 ; 0,0,0,1 ; 0,(-F/(l*M)),(((M+m)*-g)/(l*M)),0 ];
B = [0 ; (1/M) ; 0 ; (-1/(l*M))];
C = [1,0,0,0];
D = 0;

An=[A zeros(4,1) ; -C 0];
Bn=[B ; 0];
Cn=[1 0];

Ao=A';
Bo=C';
Co=B';
%comprobamos CONTROLABILIDAD del sistema
M_controlable=[B A*B A*A*B A*A*A*B];
rank(M_controlable)% nos da 4 por lo tanto es controlable


%CALCULO DEL CONTROLADOR K
%para el calculo del mismo se utiliza el metodo LQR para lo cual definimos

%Q=diag([1/10 1 0.0001 100 1]); R=100; %segundo caso
Q=diag([10000 1000 0.1 10 1]);    R=10000;
Q=1*diag([0.01   1   0.61    0.0086    0.0003]);    R=100; %por iteracion
Kn=lqr(An,Bn,Q,R);
K=Kn(1:4)
Ki=Kn(5)
eig(A-K*B); 

%Qo=diag([1 1 1 1]); R=1; %segundo caso
Qo=1e2*diag([0.2128    0.0039    5.7805    0.1200]);    R_o=0.1; %iteracion
Ko=lqr(Ao,Bo,Qo,R)

%integracion
tf=200; dt=1*10^-3; t=0:dt:(tf-dt); n=round(tf/dt);
d_i=0;           %posicion  delta   inicial
%phi_i=0;         %angulo    phi     inicial   ???
phi_i=pi;        %angulo    phi     inicial   ???
referencia1=10;  %posicion delta de referencia 1
referencia2=0;   %posicion delta de referencia 2
%generacion de funcion a usar
Ref=zeros(1,n);  masa=zeros(1,n);
for j=1:1:n-1
    if (j < (n-1)/2 )
        Ref(j)=referencia1;
        masa(j)=m;
    else (j >= (n-1)/2 )
        Ref(j)=referencia2;
        masa(j)=m*10;
    end
end
% figure
% subplot(2,1,1)
% plot(t,Ref);title('referencia de posicion');
% grid on
% subplot(2,1,2)
% plot(t,masa);title('masa variable');
% grid on

%iteracion
X_op=[0 ; 0 ; pi ; 0];
X=zeros(4,n);
X_in=[d_i 0 phi_i 0]; %[delta delta_p phi phi_p] 
X(1,1)=X_in(1);   %delta    inicial
X(2,1)=X_in(2);   %delta_p  inicial
X(3,1)=X_in(3);   %phi      inicial
X(4,1)=X_in(4);   %phi_p    inicial
psi(1)=0;         %psi      inicial

Xhat=zeros(4,n);
Xhat_in=[0 0 0 0]; %[deltahat deltahat_p phihat phihat_p] 
Xhat(1,1)=Xhat_in(1);   %deltahat    inicial
Xhat(2,1)=Xhat_in(2);   %deltahat_p  inicial
Xhat(3,1)=Xhat_in(3);   %phihat      inicial
Xhat(4,1)=Xhat_in(4);   %phihat_p    inicial

Xc=zeros(4,n);
Xc_in=[d_i 0 phi_i 0]; %[delta delta_p phi phi_p] 
Xc(1,1)=X_in(1);   %delta    inicial
Xc(2,1)=X_in(2);   %delta_p  inicial
Xc(3,1)=X_in(3);   %phi      inicial
Xc(4,1)=X_in(4);   %phi_p    inicial
psic(1)=0;         %psi      inicial

U(1)=0;
Ucom(1)=0;
% posicion de referencia a la cual se pretende llegar 

for i=1:1:n-1
    X_a=X(:,i);     %[delta ; delta_p ; phi ; phi_p ]
    Xhat_a=X(:,i);  %[deltahat ; deltahat_p ; phihat ; phihat_p ]
    Y=C*X_a;
    Yhat=C*Xhat_a;
    err=Y-Yhat;
    psi_p=Ref(i)-C*(X_a);
    psi(i+1)=psi(i)+dt*psi_p;
    Ua=-K(2:4)*(Xhat_a(2:4)-X_op(2:4))-K(1)*(X_a(1))-Ki*psi(i+1);
    Ua=-K*(Xhat_a-X_op)-Ki*psi(i+1);
    U=[U Ua];

    Xp_a1=X_a(2);
    Xp_a2=-F*X_a(2)/M-masa(i)*g*(X_a(3)-pi)/M+Ua/M;
    Xp_a3=X_a(4);
    Xp_a4=-F*X_a(2)/(M*l)-g*(M+masa(i))*(X_a(3)-pi)/(M*l)-Ua/(M*l);
    Xp_a=[Xp_a1 ; Xp_a2 ; Xp_a3 ; Xp_a4];
    Xf= X_a+ dt*(Xp_a); % Realizamos la integracion de euler y actualizamos matriz X
    X(:,i+1)=Xf;
    
    %OBSERVADOR
    Xhatp_a=err*Ko'+A*(X_a-X_op)+Ua*B;
    Xhatf=Xhat_a+dt*Xhatp_a;
    Xhat(:,i+1)= Xhatf;
    
    %COMPARACION
    Xc_a=X(:,i);    %[deltac ; deltac_p ; phic ; phic_p ]
    Yc=C*Xc_a;
    psic_p=Ref(i)-Yc;
    psic(i+1)=psic(i)+dt*psic_p;
    Uc=-K*Xc_a+Ki*psic(i+1);
    Ucom=[Ucom Uc];
    
    Xcp_a1=Xc_a(2);
    Xcp_a2=-F*Xc_a(2)/M-masa(i)*g*(Xc_a(3)-pi)/M+Uc/M;
    Xcp_a3=Xc_a(4);
    Xcp_a4=-F*Xc_a(2)/(M*l)-g*(M+masa(i))*(Xc_a(3)-pi)/(M*l)-Uc/(M*l);
    Xcp_a=[Xcp_a1 ; Xcp_a2 ; Xcp_a3 ; Xcp_a4];
    Xcf= Xc_a+ dt*(Xcp_a); % Realizamos la integracion de euler y actualizamos matriz X
    Xc(:,i+1)=Xcf;
    
    
    
end
figure
subplot(2,1,1);
plot(t,Ref,'k');
grid on
hold on
plot(t,Xc(1,:),'g');
plot(t,X(1,:),'r');title('distancia delta con d_i=0 y d_ref=10');xlabel('tiempo[s]');ylabel('posicion[m]');

legend('Referencia','Sin Observador','Con Observador');
% subplot(2,2,3);
% plot(t,U,'b');title('accion de control');


subplot(2,1,2);
%plot(t,Ref,'k');
grid on
hold on
plot(t,Xc(3,:),'g');
plot(t,X(3,:),'r');title('angulo phi con phi_i=0.03 ');xlabel('tiempo[s]');ylabel('angulo[rad]');
legend('Sin Observador','Con Observador');
% figure
% 
% plot(t,U,'b');title('accion de control');
