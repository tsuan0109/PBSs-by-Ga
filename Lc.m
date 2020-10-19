load('FitDispersive.mat');
g = 0.78;
KTM = KappaTM(g,1.55);
KTE = KappaTE(g,1.55);
LTM = pi/(2*KTM);
LTE = pi/(2*KTE);