clc; clear ; close all
%{
---------------------------------------------------------------------------
Considerar que no puede medirse la corriente y sólo pueda medirse el ángulo, por lo 
que debe implementarse un observador. Obtener la simulación en las mismas condiciones que en 
el punto anterior, y superponer las gráficas para comparar.
---------------------------------------------------------------------------
%}
%este es el sistema con las 3 variables de estado calculadas con el sistema
%observado
%DEFINO PARAMETROS
LAA = 366*10^-6;
J = 5*10^-9;
RA = 55.6;
Bm = 0;
Ki = 6.49*10^-3;
Km = 6.53*10^-3;

%DEFINO MATRICES
%X=[ia ; tita ; w];
A=[-RA/LAA 0 -Km/LAA  ; 0 0 1 ; Ki/J 0 -Bm/J];
B=[1/LAA; 0; 0];
C=[0 1 0];
D=[0];
%Dejamos explicitas las matrices del sistema dual qeu es el "observado"
An=[A zeros(3,1); -C 0];
Bn=[B ; 0];
Cn=[1 0];
Ao=A';
Bo=C';
Co=B';

%CALCULO DEL CONTROLADOR K
%para el calculo del mismo se utiliza el metodo LQR para lo cual definimos
Q=diag([1 1 1 1]); R=1;
K4=lqr(An,Bn,Q,R);
K=K4(1:3);
K_i=K4(4);
Qo=diag([1 1 1]); Ro=1;
Ko=lqr(Ao,Bo,Qo,Ro);


%implementacion de funciones a usar
tf=5; dt=1*10^-5; t=0:dt:(tf-dt); periodo=0.6;%[seg]
torq=1.15*10^-3;

Ref=pi/2*square(2*pi*t/periodo);%funcion de referencia que varia entre pi/2 y -pi/2
TL=torq/2*square(2*pi*t/periodo)+torq/2;%Funcion torque que varia entre 0 y 1.15*10^-3

%iteracion
n=round(tf/dt);
X=zeros(3,n);
X(1,1)=0; %ia inicial
X(2,1)=0; %tita inicial
X(3,1)=0; %w inicial
psi(1)=0; %psi inicial
Xhat=zeros(3,n);
Xhat(1,1)=0; %ia_hat inicial
Xhat(2,1)=0; %tita_hat inicial
Xhat(3,1)=0; %wr_hat inicial

%DEFINO CONDICIONES INICIALES

for i=1:1:n-1
    X_a=[X(1,i); X(2,i) ; X(3,i)];%[ia ; tita ; w]
    Xhat_a=[Xhat(1,i) ; Xhat(2,i) ; Xhat(3,i)];%[ia_hat ; tita_hat ; w_hat]
    Y=C*X_a;
    Yhat= Co*Xhat_a;
    psi_p=Ref(i)-Y;
    psi(i+1)=psi(i)+psi_p*dt;
    U=-K(2:3)*X_a(2:3)+K(1)*Xhat_a(1)+K_i*psi(i+1);% U estimando parametros
    U=-K*Xhat_a+K_i*psi(i+1); % U con las 3 variables de estado estimadas
    
    
    Xp_1=-RA/LAA*X_a(1)-Km/LAA*X_a(3)+1/LAA*U;  %ia_p
    Xp_2= X_a(3);                               %tita_p
    Xp_3=Ki/J*X_a(1)-Bm/J*X_a(3)-1/J*TL(i);%    %W_p
    Xp_a=[Xp_1 ; Xp_2 ; Xp_3];
    % Realizamos la integracion de euler y actualizamos matriz X
    Xf= X_a+ dt*Xp_a; 
    
    X(1,i+1)=Xf(1);
    X(2,i+1)=Xf(2);
    X(3,i+1)=Xf(3);
    %Observador
    err=(Y-Yhat);
    Xhat_p=U*Bo+Ko*err+Ao*Xhat_a;
    Xhatf=Xhat_a + dt*Xhat_p;
    
    Xhat(1,i+1)=Xhatf(1);
    Xhat(2,i+1)=Xhatf(2);
    Xhat(3,i+1)=Xhatf(3);
    
end
%ploteo de entrada con ganancia de prealimentacion U  y perturbacion TL
% figure
% subplot(2,1,1)
% hold on; grid on;
% plot(Ref);title('Referencia de Entrada');xlabel('tiempo[s]');ylabel('angulo');
% subplot(2,1,2)
% hold on;grid on;
% plot(TL);title('Torque de perturbación');xlabel('Tiempo');ylabel('Torque');


figure
plot(t,psi);
% 
% figure
% plot(t,Ref);
% grid on
% hold on
% 
% plot(t,X(2,:));title('angulo tita con observador e integrador');xlabel('tiempo[s]');ylabel('angulo[rad]');legend('REF','var de estado')







