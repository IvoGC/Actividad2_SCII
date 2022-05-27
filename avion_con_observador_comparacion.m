clc; clear ; close all
%{
---------------------------------------------------------------------------
Asumiendo que sólo puede medirse la altura, calcular un observador con los polos 
ubicados en ui=-50;-30;-1+-j. y realizar un gráfico comparativo que permita analizar el 
comportamiento de las variables de estado reales y las observadas, a lo largo del tiempo de 
simulación en cada caso, superpuestas en una sola gráfica.
---------------------------------------------------------------------------
%}
%DEFINO PARAMETROS
w=3; a=0.05; b=5; c=100;

%DEFINO MATRICES
%X=[alfa phi phi_p h]
A = [ -a,a,0,0 ; 0,0,1,0 ; w^2,-w^2,0,0 ; c,0,0,0 ];
B = [ 0; 0 ; w^2*b ; 0 ];
C = [ 0,0,0,1 ];
D = 0;

Ao=A';
Bo=C';
Co=B';


%CALCULO DEL CONTROLADOR K
%en este caso tenemos el valor de los polos deseados y en base a esto
%obtenemos el valor del controlador por lo tanto utilizando la formula de
%akerman lo impementamos
p1=-15+15i;
p2=-15-15i;
p3=-0.5+0.5i;
p4=-0.5-0.5i;

%syms p1 p2 p3 p4 s
%expand((s-p1)*(s-p2)*(s-p3)*(s-p4))
a0=1;
a1=-p1-p2-p3-p4;
a2=p1*p2+p1*p3+p1*p4+p2*p3+p2*p4+p3*p4;
a3=-p1*p2*p3-p1*p2*p4-p1*p3*p4-p2*p3*p4;
a4=p1*p2*p3*p4;
phi_A=a0*A^4+a1*A^3+a2*A^2+a3*A^1+a4*A^0;

AUX=[B A*B A*A*B A*A*A*B];
K=[0 0 0 1]*inv(AUX)*phi_A;



%ahora para el observador
po1=-50;
po2=-30;
po3=-1+i;
po4=-1-i;

polos_o=[po1 po2 po3 po4];
Ko=place(Ao, Bo, polos_o);

%Ganancia de prealimentacion
G=-inv(C*inv(A-B*K)*B);

%simulacion Para altura inicial 
tf=70; dt=1*10^-3; t=0:dt:(tf-dt); 
h_o=500;       %Altura inicial
referencia=-100; %Valor de establecimiento
n=round(tf/dt);

X=zeros(4,n);
X_in=[0 0 0 h_o]; %[alfa phi phi_p h] ho=4000
X(1,1)=X_in(1);   %alfa   inicial
X(2,1)=X_in(2);   %phi    inicial
X(3,1)=X_in(3);   %phi_p  inicial
X(4,1)=X_in(4);   %h      inicial


Xhat=zeros(4,n);
Xhat_in=[0 0 0 0]; %[alfahat phihat phi_phat hhat]
Xhat(1,1)=Xhat_in(1);   %alfahat   inicial
Xhat(2,1)=Xhat_in(2);   %phihat    inicial
Xhat(3,1)=Xhat_in(3);   %phihat_p  inicial
Xhat(4,1)=Xhat_in(4);   %hhat      inicial


Xc_in=[0 0 0 h_o]; %[alfac phic phic_p hc] ho=4000
Xc(1,1)=Xc_in(1);   %alfac   inicial
Xc(2,1)=Xc_in(2);   %phic    inicial
Xc(3,1)=Xc_in(3);   %phic_p  inicial
Xc(4,1)=Xc_in(4);   %hc      inicial

Ua(1)=0;
Uc(1)=0;
%DEFINO CONDICIONES INICIALES
Ref=referencia*ones(1,n);% altura de referencia a la cual se pretende llegar (altura referenciada al ho)
for i=1:1:n-1
    X_a=X(:,i);         %[alfa ; phi ; phi_p ; h]
    Xhat_a=Xhat(:,i);   %[alfahat ; phihat ; phihat_p ; hhat]
    Y=C*X_a;
    Yhat=C*Xhat_a;
    err=Y-Yhat;
    %U=-K(1:3)*Xhat_a(1:3)-K(4)*X_a(4)+Ref(i)*G;
    %entrada U estimando otros parametros y midiendo altura
    %U=-K*Xhat_a+Ref(i)*G;
    U=-K(1:3)*Xhat_a(1:3)-K(4)*X_a(4)+Ref(i)*G;
    
    %entrada U estimando los 3 parametros
    Ua=[Ua U];
    
    Xp_1=a*(X_a(2)-X_a(1));         %alfa_p
    Xp_2= X_a(3);                   %phi_p
    Xp_3=-w^2*(X_a(2)-X_a(1)-b*U);  %phi_pp
    Xp_4=c*X_a(1);                  %h_p
    Xp_a=[Xp_1 ; Xp_2 ; Xp_3 ; Xp_4];
    Xf= X_a+ dt*Xp_a; % Realizamos la integracion de euler y actualizamos matriz X
    X(:,i+1)=Xf;
    
    %OBSERVADOR
    Xhat_p=U*B+Ko'*err+A*Xhat_a;
    Xhatf= Xhat_a+ dt*Xhat_p; % Realizamos la integracion de euler y actualizamos matriz X
    Xhat(:,i+1)=Xhatf;
    
    %SIN OBSERVADOR A COMPARAR
    Xc_a=Xc(:,i);
    Yc=C*Xc;
    Uc=-K*Xc_a+Ref(i)*G;
%    Uc=[Uc U];
    Xcp_1=a*(Xc_a(2)-Xc_a(1));         %alfac_p
    Xcp_2= Xc_a(3);                   %phic_p
    Xcp_3=-w^2*(Xc_a(2)-Xc_a(1)-b*Uc);  %phic_pp
    Xcp_4=c*Xc_a(1);                  %hc_p
    Xcp_a=[Xcp_1 ; Xcp_2 ; Xcp_3 ; Xcp_4];
    Xcf= Xc_a+ dt*Xcp_a; % Realizamos la integracion de euler y actualizamos matriz X
    Xc(:,i+1)=Xcf;
    
end
figure
subplot(2,2,1);
grid on
hold on
plot(t,Xc(4,:),'g');
plot(t,Ref,'k');
plot(t,X(4,:),'r');title('comparacion altura h');xlabel('tiempo[s]');ylabel('posicion[m]');
legend('sin observador','Ref','con observador')

subplot(2,2,2);
grid on
hold on
plot(t,Xc(1,:),'g');
%plot(t,Ref,'k');
plot(t,X(1,:),'r');title('comparacion alfa');xlabel('tiempo[s]');ylabel('angulo[RAD]');
legend('sin observador','con observador')

subplot(2,2,3);
grid on
hold on
plot(t,Xc(2,:),'g');
%plot(t,Ref,'k');
plot(t,X(2,:),'r');title('comparacion phi');xlabel('tiempo[s]');ylabel('angulo[RAD]');
legend('sin observador','con observador')

subplot(2,2,4);
grid 
hold on
plot(t,Xc(3,:),'g');
%plot(t,Ref,'k');
plot(t,X(3,:),'r');title('comparacion phi_p ');xlabel('tiempo[s]');%ylabel('vel angular[]');
legend('sin observador','con observador')






