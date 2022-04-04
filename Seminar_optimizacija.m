t = 0:0.01:80; %vrijeme simulacije
H=200; %ukupna visina
h0=2; %visina lifta 
h1=1; % razlika u visini koluta 
h2=1; %visina protuutega
ro=0.1; % masa kabela po metru
m1=60; %masa lifta bez ljudi
m2=90; %masa protuutega
R=0.25; %radijus koluta
J=0.1; %moment inercije kotača zajedno sa inercijom osovine motora + reduktor
b1=0.5; % trenje vodilica dizala
b2=0.45; %trenje vodilica protuutega;
g=9.81; %gravitacija
x10=0; % promjenjiva vrijednost početnog položaja,a maksimalno je x10max = H-h1-h0
x20=H-h1-x10; %početni položaj protuutega
dm=20; % masa ljudi
%proracun koeficijenata u diferencijalnoj jednadžbi
% a*theta''+b*theta'+c*theta=d+upravljacki signal(moment na kolutu)
a = J+R^2*(m1+m2+dm)+R^2*ro*(2*H-x10-x20-h0-h1-h2);
b = R^2*(b1+b2);
c = -2*R^2*ro*g;
d = R*g*(m2-m1-dm+ro*(x10-x20+h0-h1-h2));
Tlstac = -d; % konstantna sila koja je potrebna da lift miruje u x=0 sve derivacije i theta=0

% matrice varijabli stanja 
A = [0 1;-c/a -b/a];
B = [0;1/a];
C = [1 0];
D = 0;
sys = ss(A,B,C,D); % sustva u otvorenoj petlji
SYSd=c2d(sys,0.01,'zoh'); %diskretizacija (Eulerova aproksimacija?)
Ad=SYSd.A; %izlučivanje matrica A i B (fi i psi) za iteriranje
Bd=SYSd.B;

%%% odziv sustava s neg. povrtanom vezom na jedinicnu pobudu  
T = d+zeros(length(t),1); % referentnom signalu T dodao sam d, zbog
T(t>=0) = 1+d;            % konstante koja djeluje zbog gravitacije itd.

figure(1)
SYSdf=feedback(SYSd,1); %zatvaranje sustava sa neg. povratnom vezom
lsimplot(SYSdf,T,t)     %crtanje 
title('Odziv sustava sa povratnom vezom na jedinični step')
ylabel('Položaj lifta')
xlabel('Vrijeme')
grid on
[Y,t0,X]=lsim(SYSdf,T,t,[0;0]); %spremanje izlaznog vektora Y(x(1),položaj u radijanima),
figure(2)                      % vektora vremena t0 i matrice stanja(2 vektora stanja)
plot(X(:,1),X(:,2)) %fazna ravnina
grid on
title('Fazna ravnina,žarište sustava bez regulatora')
xlabel('Položaj')
ylabel('Brzina')

%vjerovatno postoji fleksibilniji nacin za generiranje vektora brzine,
% praktički se za svaki kat može nacrtati graf i po njemu izračunati 
% položaj i akceleracija

% % željeni brzina,akceleracija i položaj u rad/s,rad/s^2 i rad

%brzina
v=zeros(length(t),1); 
v(t<6)=1.6*t(t<6);   
v(t>24 & t<=30)=flip(v(t<6));
nagib = flip(v(t<6));
v(t>=6 & t<=24)=max(v);
v(t>30 & t<35)=0;
v(t>=35 & t<41)=-1.6*t(t<6);
v(t>=41 & t<=51)=-max(v);
v(t>51 & t<=57)=-nagib;

%akceleracija
akc=gradient(v,t);

r=zeros(length(t),1);
r0=0;
v0=0;
a0=0;

%položaj
for i=1:length(t)-1
r(i) = r0 + (v0+v(i+1)-v(i))*(t(i+1)-t(i)) + 0.5*akc(i)*(t(i+1)-t(i))^2;
r0=r(i);
v0=v(i);
end
r(end)=r(8000);

% priprema referentnog vektora stanja usporedbu kod računanja cijene
xref=[r';v'];

xk =[0;0]; % početni polozaj i brzina u radijanima,rad/s
XK =[0;0]; %isto za petlju u kojoj primjenjujemo nađene najbolje
           %upravljacke signale, pomičemo horizont i pozivamo fukciju u
           %kojoj primjenjujemo seriju upravljački signala,a za svaki
           %nađemo sva stanja u horizontu,proračunamo cijenu za taj signal
           % i po najmanjoj cijeni odabiremo najbolji za taj trenutak

stanja=zeros(2,length(t)); %alokacija memorije
CTRL=zeros(1,length(t));

for z=1:length(t) % petlja do krja simulacije
    %poziv funkcije koja vraca najbolji upravljacki signal ctrl
    [ctrl,FJ,x_DOB,x_REF0] = najbolji_upravljacki(z,xk,xref,Ad,Bd,d); 
    
    ctrlf=filtracija(z,ctrl);
     % spremanje najboljeg upravljackog signala u vektor za kasniju simulaciju
    CTRL(z)=ctrlf+d;
    
    %iteracija pocevsi od XK ,a controlni signal se mijenja u svakoj iteraciji
    %novo stanje=A*staro stanje + B*upravljacki signal
    %nisam siguran u dodavanja konstante d tu.
    stanja(:,z)=Ad*XK(:,z)+Bd*(ctrlf+d);
    
    % novo stanje je staro
    XK(:,z+1)=stanja(:,z);
    
    %spremanje novog stanja za ulaz u funkciju
    %kao pocetno iz kojega se simuliraju nove
    %trajektorije za različite upravljačke signale
    xk=stanja(:,z);
end

figure(3)
% usporedba zeljenog odziva brzine i akceleracije sa dobivenim
plot(t,v*R,'b')
hold on
plot(t,stanja(2,:)*R,'r')
plot(t,akc*R,'g')
plot(t,gradient(stanja(2,:),t)*R,'k')
grid on
title('Usporedba željenih i dobivenih brzina i akceleracija u m/s i m/s^2')
ylabel('Brzina / Akceleracija')
xlabel('Vrijeme')
hold off

% crtanje položaja 
figure(4)
plot(t,R*r,'r')
hold on
[P]=lsim(SYSd,CTRL,t);
plot(t,P*R,'b')
plot(t,CTRL,'k')

%dodavanje izglađenog kontrolnog signala i crtanje odziva sustava na takav
%signal

%realniji_CTRL = smoothdata(CTRL,'gaussian',20); %izglađivanja zbog realnije prikaza
%plot(t,realniji_CTRL,'Linewidth',2)
%[Po]=lsim(SYSd,realniji_CTRL,t);
%plot(t,Po*R)
grid on
title('Usporedba željene sa dobivenom putanjom lifta i upravljački signal')
ylabel('Položaj')
xlabel('Vrijeme')
% crtanje dijela željene i dobivene trajektorije za neki upravljački
% signal u 5000-toj iteraciji u simulaciji i prikaz proračuna cijene za taj
% upravljački signal

figure(5)
plot(x_DOB(1,:),x_DOB(2,:),'*b')
hold on
plot(x_REF0(1,:),x_REF0(2,:),'xr')
title('Usporedba dijela željene i dobivene trajektorije' )
xlabel('Položaj')
ylabel('Brzina')
text(114.1,-7,'Cijena je druga norma, odnosno '+string(FJ));
grid on
plot([x_REF0(1,:);x_DOB(1,:)],[x_REF0(2,:);x_DOB(2,:)])

%%% za simulink 
ts=timeseries(CTRL,t);
