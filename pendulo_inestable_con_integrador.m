clc; clear ; close all
%{
---------------------------------------------------------------------------
Calcular un controlador que haga evolucionar al péndulo en el equilibrio inestable, 
partiendo de una condición inicial nula en el desplazamiento y termine en 10 metros manteniendo la 
vertical. Determinar el ángulo máximo que puede alejarse de la vertical en t=0 para que el sistema 
cumpla el objetivo de control. 
---------------------------------------------------------------------------
%}
%en este caso con integrador
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
%Matrices extendidas para el integrador
An=[A zeros(4,1); -C 0];
Bn=[B ; 0];
Cn=[1 0];


%CALCULO DEL CONTROLADOR K
%para el calculo del mismo se utiliza el metodo LQR para lo cual definimos
Q=diag([1/50 1/8000 1/500 1/1000 1]); R=0.01;
%Q=diag([1 1 1 1]); R=1;
K5=lqr(An,Bn,Q,R);
K=K5(1:4)
Ki=-K5(5)
eig(A-K*B); 



%integracion

tf=10; dt=1*10^-3; t=0:dt:(tf-dt); 
d_i=0;          %posicion delta inicial
phi_i=0.03;        %angulo   phi   inicial
referencia=10;  %posicion delta de referencia
%iteracion
n=round(tf/dt);
X=zeros(4,n);
X_in=[d_i 0 phi_i 0]; %[delta delta_p phi phi_p] 
X(1,1)=X_in(1);   %delta    inicial
X(2,1)=X_in(2);   %delta_p  inicial
X(3,1)=X_in(3);   %phi      inicial
X(4,1)=X_in(4);   %phi_p    inicial
psi(1)=0;         %psi      inicial

U(1)=0;
Ref=referencia*ones(1,n);% posicion de referencia a la cual se pretende llegar 

for i=1:1:n-1
    X_a=X(:,i);%[delta ; delta_p ; phi ; phi_p ]
    psi_p=Ref(i)-C*X_a;
    psi(i+1)=psi(i)+dt*psi_p;
    Ua=-K*X_a+Ki*psi(i+1);
    U=[U Ua];
    
    Xp_a=A*X_a+B*Ua;
    Xf= X_a+ dt*Xp_a; % Realizamos la integracion de euler y actualizamos matriz X
    X(:,i+1)=Xf;
end
figure
subplot(2,1,1);
plot(t,Ref,'k');
grid on
hold on
plot(t,X(1,:),'r');title('distancia delta con d_i=0 y d_ref=10');xlabel('tiempo[s]');ylabel('posicion[m]');
legend('Referencia','posicion')
% subplot(2,2,3);
% plot(t,U,'b');title('accion de control');


subplot(2,1,2);
%plot(t,Ref,'k');
grid on
hold on
plot(t,X(3,:),'r');title('angulo phi con phi_i=0.03 ');xlabel('tiempo[s]');ylabel('angulo[rad]');
legend('angulo')
figure

plot(t,U,'b');title('accion de control');

