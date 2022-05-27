clc; clear ; close all
%{
---------------------------------------------------------------------------
Incorporar un observador para el caso en que sólo pueda medirse el 
desplazamiento delta, repetir las simulaciones para las condiciones anteriores 
y graficar los resultados en gráficas superpuestas.
---------------------------------------------------------------------------
%}
%En este caso vamos usar el caso sin integrador ya que tiene mejor
%respuesta que el caso con el integrador y no es necesario el mismo. 
%vamos a generar el sistema con observador y un sistema de comparacion para
%poder plotearlo en simultaneo
%DEFINO PARAMETROS
m = 0.1; F = 0.1; l = 0.6; g = 9.8; M = 0.5; 

%DEFINO MATRICES
%X=[delta delta_p phi phi_p]
%Sistema linealizado en el punto Xo
X_o=[0 0 0 0];%equilibio inestable phi=0
A = [0,1,0,0 ; 0,(-F/M),(-m*g/M),0 ; 0,0,0,1 ; 0,(F/(l*M)),(((M+m)*g)/(l*M)),0 ];
B = [0 ; (1/M) ; 0 ; (-1/(l*M))];
C = [1,0,0,0];
D = 0;
Ao=A';
Bo=C';
Co=B';

%CALCULO DEL CONTROLADOR K
%para el calculo del mismo se utiliza el metodo LQR para lo cual definimos
Q=diag([1/300 1/80000 1/5000 1/100]); R=0.01;
%Q=diag([1 1 1 1]); R=1;
K=lqr(A,B,Q,R)
eig(A-K*B); 

Qo=diag([100 1/100 1 1/500]); Ro=0.01;
Ko=lqr(Ao,Bo,Qo,Ro)

%Ganancia de prealimentacion
G=-inv(C*inv(A-B*K)*B);

%integracion
tf=30; dt=1*10^-3; t=0:dt:(tf-dt); 
d_i=0;              %posicion delta inicial
phi_i=-0.03         %angulo   phi   inicial
referencia=10;      %posicion delta de referencia
%iteracion
n=round(tf/dt);
X=zeros(4,n);
X_in=[d_i 0 phi_i 0]; %[delta delta_p phi phi_p] 
X(1,1)=X_in(1);   %delta    inicial
X(2,1)=X_in(2);   %delta_p  inicial
X(3,1)=X_in(3);   %phi      inicial
X(4,1)=X_in(4);   %phi_p    inicial

Xhat=zeros(4,n);
Xhat_in=[0 0 0 0]; %[deltahat deltahat_p phihat phihat_p] 
Xhat(1,1)=Xhat_in(1);   %deltahat       inicial
Xhat(2,1)=Xhat_in(2);   %deltahat_p     inicial
Xhat(3,1)=Xhat_in(3);   %phihat         inicial
Xhat(4,1)=Xhat_in(4);   %phihat_p       inicial

Xc=zeros(4,n);
Xc_in=[d_i 0 0 phi_i]; %[deltac deltac_p phic phic_p] 
Xc(1,1)=Xc_in(1);   %deltac     inicial
Xc(2,1)=Xc_in(2);   %deltac_p   inicial
Xc(3,1)=Xc_in(3);   %phic       inicial
Xc(4,1)=Xc_in(4);   %phic_p     inicial


U(1)=0;
Uc(1)=0;
Ref=referencia*ones(1,n);% posicion de referencia a la cual se pretende llegar 

for i=1:1:n-1
    X_a=X(:,i);         %[delta ; delta_p ; phi ; phi_p ]
    Xhat_a=Xhat(:,i);   %[deltahat deltahat_p phihat phihat_p] 
    Y=X_a*C;
    Yhat=Xhat_a*C;
    err=Y-Yhat;
    %Ua=-K*X_a+Ref(i)*G;             %U del sistema normal
    %Ua=-K*Xhat_a+Ref(i)*G;          %U del sistema observado
    Ua=-K(2:4)*Xhat_a(2:4)-K(1)*X_a(1)+Ref(i)*G;%U para el sistema mixto
    U=[U Ua];
    Xp_a=A*X_a+B*Ua;
    X(:,i+1)= X_a+ dt*Xp_a;
    % %Xf= X_a+ dt*Xp_a; % Realizamos la integracion de euler y actualizamos matriz X
    % %X(:,i+1)=Xf;
    %OBSERVADOR
    Xhat_p=err*Ko'+A*Xhat_a+Ua*B;
    Xhat(:,i+1)=Xhat_a+dt*Xhat_p;
    Xhatf=Xhat_a+dt*Xhat_p;
    Xhat(:,i+1)=Xhatf;
    
    %sistema de COMPARACION sin observador
    Xc_a=Xc(:,i);       %[deltac ; deltac_p ; phic ; phic_p ]
    Uc_a=-K*Xc_a+Ref(i)*G;
    Uc=[Uc Uc_a];
    Xc_p=A*Xc_a+B*Uc_a;
    Xcf=Xc_a+ dt*Xc_p;
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
