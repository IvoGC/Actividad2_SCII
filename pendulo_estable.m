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

%comprobamos CONTROLABILIDAD del sistema
M_controlable=[B A*B A*A*B A*A*A*B];
rank(M_controlable)% nos da 4 por lo tanto es controlable


%CALCULO DEL CONTROLADOR K
%para el calculo del mismo se utiliza el metodo LQR para lo cual definimos

%Q=diag([1/300 1/80000 1/5000 1/100]); R=0.01; %segundo caso
%Q=diag([0.1 1 1/800000 1/50]); R=10;
Q=1*diag([1000 1 1000 10]);    R=.1;
K=lqr(A,B,Q,R);
eig(A-K*B) %

%Ganancia de prealimentacion
G=-inv(C*inv(A-B*K)*B);

%integracion
tf=100; dt=1*10^-3; t=0:dt:(tf-dt); n=round(tf/dt);
d_i=0;           %posicion delta inicial
phi_i=pi;        %angulo   phi   inicial
referencia1=10;  %posicion delta de referencia 1
referencia2=0;   %posicion delta de referencia 2
%generacion de funcion a usar
Ref=zeros(1,n);  masa=zeros(1,n);
for j=1:1:n-1
    if (j < (n-1)/2 )
        Ref(j)=10;
        masa(j)=m;
    else (j >= (n-1)/2 )
        Ref(j)=0;
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

U(1)=0;
% posicion de referencia a la cual se pretende llegar 

for i=1:1:n-1
    X_a=X(:,i);%[delta ; delta_p ; phi ; phi_p ]
    Ua=-K*(X_a-X_op)+Ref(i)*G;
    U=[U Ua];
%     delta_pp = (-F*y(2,i)/M)-(m*g*(y(3,i)-pi)/M)+(u/M);
%     phi_pp   = (-F*y(2,i)/(l*M))-(g*(M+m)*(y(3,i)-pi)/(M*l))-u/(l*M);
    Xp_a1=X_a(2);
    Xp_a2=-F*X_a(2)/M-masa(i)*g*X_a(3)/M+Ua/M;
    Xp_a3=X_a(4);
    Xp_a4=-F*X_a(2)/(M*l)-g*(M+masa(i))*X_a(3)/(M*l)-Ua/(M*l);
    
    Xp_a=[Xp_a1 ; Xp_a2 ; Xp_a3 ; Xp_a4];
    Xf= X_a+ dt*(Xp_a-X_op); % Realizamos la integracion de euler y actualizamos matriz X
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
% figure
% 
% plot(t,U,'b');title('accion de control');

