function A = CMT(p)
[c,ceq] = Constraint(p);
load('FitDispersive.mat');
if c(1)>0 || c(2)>0 || c(3)>0 || c(4)>0
    A = 100;
else

%網格數
L = 48;   %這邊L是指扣掉頭尾的長度，實際長度為48
dL = 1;
dl = 0.1;  
N0 = L/dL+1;                   
N1 = 1+L/2/dL;           %中間那格,其中中間設定為0.48 0.48                   
N = L/dl+1;

%p = [1.41315001286629	1.32630002573257	1.23945003859886	1.15260005146514	1.08238189563708	1.01216373980903	0.941945583980969	0.871727428152912	0.785516662589317	0.699305897025721	0.613095131462126	0.526884365898530	0.479177422987134	0.431470480075739	0.383763537164343	0.336056594252947	0.334665700566894	0.333274806880841	0.331883913194787	0.330493019508734	0.312896058703437	0.295299097898139	0.277702137092842	0.260105176287544 0.483959590001032	0.487919180002063	0.491878770003095	0.495838360004127	0.502576523409380	0.509314686814633	0.516052850219886	0.522791013625139	0.512381259841181	0.501971506057224	0.491561752273266	0.481151998489308	0.485001741203411	0.488851483917514	0.492701226631616	0.496550969345719	0.506515517609466	0.516480065873213	0.526444614136960	0.536409162400707	0.522306871800530	0.508204581200354	0.494102290600177];

%Gap
% p = [1.25,1,0.75,0.45,0.45,0.45,0.48,0.48,0.48,0.48,0.48,0.48]

b = p(7:11);

w1(1:5) = linspace(0.48,b(1),5);
w1(5:9) = linspace(b(1),b(2),5);
w1(9:13) = linspace(b(2),b(3),5);
w1(13:17) = linspace(b(3),b(4),5);
w1(17:21) = linspace(b(4),b(5),5);
w1(21:25) = linspace(b(5),0.48,5);
w2(1:N1) = 0.48*2-w1(1:N1);
w1(N1+1:N0) = w2(N1-1:-1:1);
w2(N1+1:N0) = w1(N1-1:-1:1);

V1 = 0:dL:L;
V2 = 0:dl:L;
w3 = interp1(V1,w1,V2);
w4 = interp1(V1,w2,V2);        

a = p(1:6);

g(1:5) = linspace(1.5,a(1),5);
g(5:9) = linspace(a(1),a(2),5);
g(9:13) = linspace(a(2),a(3),5);
g(13:17) = linspace(a(3),a(4),5);
g(17:21) = linspace(a(4),a(5),5);
g(21:25) = linspace(a(5),a(6),5);

g1(1:N1) = g(1:N1);
g1(N1+1:N0) = g1(N1-1:-1:1);
g2 = interp1(V1,g1,V2); 


ETE = zeros(1,5);
     for k = 1:5
        %1.5:1.6 um
        x = 1.5+(k-1)*0.025;                           %從lamda = 1.5 ~.1.6
        KTE = KappaTE(g2,x*ones(1,N));
        betaTE1 = BetaTE(w3,x*ones(1,N));
        betaTE2 = BetaTE(w4,x*ones(1,N));
        dBTE = (betaTE1-betaTE2)/2;
        PTE = sqrt(KTE.^2+dBTE.^2);
        %transfer matrices
        TE11 = cos(PTE.*dl)+1j.*dBTE.*sin(PTE.*dl)./PTE;            %T11 T22 T12 T21 分別都是1*501個矩陣,因為db.K為矩陣 為什麼少一個負號???????
        TE22 = cos(PTE.*dl)-1j.*dBTE.*sin(PTE.*dl)./PTE;
        TE12 = 1j.*KTE.*sin(PTE.*dl)./PTE;
        TE21 = 1j.*KTE.*sin(PTE.*dl)./PTE;
        % Put Matrices in a cell
        TTE = cell(1,N);                             %T裡面有501個Cell,每個cell裡面包含4個元素
        
        for i = 1:N
                TTE{i} = zeros(2,2);
                TTE{i}(1,1) = TE11(1,i);
                TTE{i}(2,2) = TE22(1,i);
                TTE{i}(1,2) = TE12(1,i);
                TTE{i}(2,1) = TE21(1,i);
        end
        UTE = [1;0];      %一個波導的input振幅為1 另一個為0
        for i = 1:N
            UTE = TTE{i}*UTE;  
        end
        ETE(1,k)=UTE(1).*conj(UTE(1));    %把11個波長算出來的output丟到E
        ETE(2,k)=UTE(2).*conj(UTE(2));    %把11個波長算出來的output丟到E
     end  
     
ETM = zeros(1,5);
     for k = 1:5
        %1.5:1.6 um
        x = 1.5+(k-1)*0.025;                           %從lamda = 1.5 ~.1.6
        KTM = KappaTM(g2,x*ones(1,N));
        betaTM1 = BetaTM(w3,x*ones(1,N));
        betaTM2 = BetaTM(w4,x*ones(1,N));
        dBTM = (betaTM1-betaTM2)/2;
        PTM = sqrt(KTM.^2+dBTM.^2);
        %transfer matrices
        TM11 = cos(PTM.*dl)+1j.*dBTM.*sin(PTM.*dl)./PTM;            %T11 T22 T12 T21 分別都是1*501個矩陣,因為db.K為矩陣 為什麼少一個負號???????
        TM22 = cos(PTM.*dl)-1j.*dBTM.*sin(PTM.*dl)./PTM;
        TM12 = 1j.*KTM.*sin(PTM.*dl)./PTM;
        TM21 = 1j.*KTM.*sin(PTM.*dl)./PTM;
        % Put Matrices in a cell
        TTM = cell(1,N);                             %T裡面有501個Cell,每個cell裡面包含4個元素
        
        for i = 1:N
                TTM{i} = zeros(2,2);
                TTM{i}(1,1) = TM11(1,i);
                TTM{i}(2,2) = TM22(1,i);
                TTM{i}(1,2) = TM12(1,i);
                TTM{i}(2,1) = TM21(1,i);
        end
        UTM = [1;0];      %一個波導的input振幅為1 另一個為0
        for i = 1:N
            UTM = TTM{i}*UTM;  
        end
        ETM(1,k)=UTM(1).*conj(UTM(1));   
        ETM(2,k)=UTM(2).*conj(UTM(2));    
     end
         e(1:5) = ETM(1,1:5);
         e(6:10) = ETE(2,1:5);        
         e = e.^2;
         A = sqrt(sum(e)/10);
end
 end
