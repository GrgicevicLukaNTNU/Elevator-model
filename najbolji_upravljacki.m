
% funkcija daje najbolji upravljački signal iz podataka o stanju iz kojega 
% počimaju iteracije xK ili xk, početku horizonta pol_horiz ili z, te
% referentu matricu po kojoj gledamo koliko odstupamo od željene brzine i
% položaja u rad i rad/s.
%također potrebno je proslijediti matrice sustava i varijablu d za dodati
%na ulaz.

function [Upr,FJ,x_DOB,x_REF0] = najbolji_upravljacki(pol_horiz,xK,xref,Ad,Bd,d)
horizont=20; %duljina horizonta = 20*0.01= 0.2 s
x_pocetni=xK; % prebacivanje vektora stanja iz kojega pocimamo u x_pocetni radi
              % reinicijalizacije nakon iteriranja do kraja horizonta
%alokacija memorije za cijenu za jedan upravljacki signal do kraja
%horizonta 
fJ=zeros(horizont,1);
%alokacija memorije za kumulativnu cijenu za sve upravljacke signale
%svaki redak je suma cijene za jedan upravljacki signal
fJk=zeros(200,1);
%alokacija memorije za primjenjene upravljačke signale, i tog vektora se
%kasnije bira najbolji,koji daje namanju razliku u stanjima
U0=zeros(200,1);

k=1;
kk=pol_horiz-1; % sluzi za pomicanje u određeni stupac u vektoru referentnih
                %stanja kako napredujemo kroz simulaciju (pol_horiz=z)
                
for U=-200:2:200 %paleta primjenjenih upravljačkih signala
U0(k)=U; % spremanje u vektor iz kojega biramo najbolji
xkk = zeros(2,horizont); %alokacija memorije za nova stanja,odnosno
                         %reinicijalizacija novih stanja u nulu nakon
                         %primjene jednom upravljackog signala i iteriranja
                         %do kraja horizonta.
for j=1:horizont   %iteriranje 
xkk(:,j) = Ad*xK(:,j) + Bd*(U+d); %nisam siguran u dodavnaje konstane d tu.
xK(:,j+1)=xkk(:,j); %novo stanje je staro stanje

if (kk+j) > 8001 % provjera ako horizont 'gleda' preko kraja simulacije 
    kk=kk-1;     %ako da,onda samo smanjimo kk ,koji je uvjetovan brojem z
end              % koji nam govori gdje smo u simulaciji.

fJ(j)=norm((xkk(:,j)-xref(:,j+kk)),2); % proračun cijene kao Euklidove norme
% odnosno prostornom udaljenosti između referentnog i vektora stanja kojega smo
% dobili iteriranjem.
end

%spremanje vektora xkk i xref za crtanje i proračun cijene
%za random upravljački signal
global x_DOB;
global x_REF0;
global FJ;
if pol_horiz==5000
x_DOB=xkk;
x_REF0=xref(:,pol_horiz:pol_horiz+19);
FJ=norm(x_DOB-x_REF0,2);
end

fJk(k,:)=sum(fJ); % sumiranje cijene za svaku granu i spremenje u k redaka
fJ=zeros(horizont,1); % reinicijalizacija vektora cijena na nulu
k=k+1;
xK=x_pocetni; %reinicijalizacija vektora starog stanja na stanje
%iz ulaza u funkciju xK(x_pocetni) za novi upravljacki signal
end
minimum=find(fJk==min(fJk)); % traženje minimuma u vektoru sumiranih cijena
%tako da formiramo vektor sa nulama koji ima 1 samo gdje je minimum

[red,stup] = size(U0(minimum)); %izbjegavanje greške ako postoje više 
                                % upravljačkih signala sa istom cijenom
if (red > 1) || (stup > 1)      %pa je izlaz iz funkcije vektor a ne skalar
    scalar = 0;
else
    scalar = 1;
end
if scalar==true
Upr = U0(minimum); % najbolji upravljački signal je na indexu gdje minimum ima 1,
%u vektoru U0 gdje su spremljeni primjenjeni upravljacki signali
else
Upr=0;
end

end

