#ifndef QUADRATURES_H
#define QUADRATURES_H

#define GET_LAG_ROOTS(d, xi, wi) switch(d) { case 5: xi = lag_p005; wi = lag_w005; break; case 10: xi = lag_p010; wi = lag_w010; break; case 15: xi = lag_p015; wi = lag_w015; break; case 20: xi = lag_p020; wi = lag_w020; break; case 25: xi = lag_p025; wi = lag_w025; break; case 30: xi = lag_p030; wi = lag_w030; break; case 35: xi = lag_p035; wi = lag_w035; break; case 40: xi = lag_p040; wi = lag_w040; break; case 45: xi = lag_p045; wi = lag_w045; break; case 50: xi = lag_p050; wi = lag_w050; break; case 55: xi = lag_p055; wi = lag_w055; break; case 60: xi = lag_p060; wi = lag_w060; break; case 65: xi = lag_p065; wi = lag_w065; break; case 70: xi = lag_p070; wi = lag_w070; break; case 75: xi = lag_p075; wi = lag_w075; break; case 80: xi = lag_p080; wi = lag_w080; break; case 85: xi = lag_p085; wi = lag_w085; break; case 90: xi = lag_p090; wi = lag_w090; break; case 95: xi = lag_p095; wi = lag_w095; break; case 100: xi = lag_p100; wi = lag_w100; break; case 105: xi = lag_p105; wi = lag_w105; break; case 110: xi = lag_p110; wi = lag_w110; break; case 115: xi = lag_p115; wi = lag_w115; break; case 120: xi = lag_p120; wi = lag_w120; break; case 125: xi = lag_p125; wi = lag_w125; break; case 130: xi = lag_p130; wi = lag_w130; break; case 135: xi = lag_p135; wi = lag_w135; break; case 140: xi = lag_p140; wi = lag_w140; break; case 145: xi = lag_p145; wi = lag_w145; break; case 150: xi = lag_p150; wi = lag_w150; break; case 155: xi = lag_p155; wi = lag_w155; break; case 160: xi = lag_p160; wi = lag_w160; break; case 165: xi = lag_p165; wi = lag_w165; break; case 170: xi = lag_p170; wi = lag_w170; break; case 175: xi = lag_p175; wi = lag_w175; break; case 180: xi = lag_p180; wi = lag_w180; break; case 185: xi = lag_p185; wi = lag_w185; break; case 190: xi = lag_p190; wi = lag_w190; break; case 195: xi = lag_p195; wi = lag_w195; break; case 200: xi = lag_p200; wi = lag_w200; break; case 205: xi = lag_p205; wi = lag_w205; break; case 210: xi = lag_p210; wi = lag_w210; break; case 215: xi = lag_p215; wi = lag_w215; break; case 220: xi = lag_p220; wi = lag_w220; break; case 225: xi = lag_p225; wi = lag_w225; break; case 230: xi = lag_p230; wi = lag_w230; break; case 235: xi = lag_p235; wi = lag_w235; break; case 240: xi = lag_p240; wi = lag_w240; break; case 245: xi = lag_p245; wi = lag_w245; break; case 250: xi = lag_p250; wi = lag_w250; break; case 255: xi = lag_p255; wi = lag_w255; break; case 260: xi = lag_p260; wi = lag_w260; break; case 265: xi = lag_p265; wi = lag_w265; break; case 270: xi = lag_p270; wi = lag_w270; break; case 275: xi = lag_p275; wi = lag_w275; break; case 280: xi = lag_p280; wi = lag_w280; break; case 285: xi = lag_p285; wi = lag_w285; break; case 290: xi = lag_p290; wi = lag_w290; break; case 295: xi = lag_p295; wi = lag_w295; break; case 300: xi = lag_p300; wi = lag_w300; break; default: xi = arma::vec(); wi = arma::vec(); }
// Gauss-Laguerre quadrature for d=5
extern arma::vec lag_p005;
extern arma::vec lag_w005;
// Gauss-Laguerre quadrature for d=10
extern arma::vec lag_p010;
extern arma::vec lag_w010;
// Gauss-Laguerre quadrature for d=15
extern arma::vec lag_p015;
extern arma::vec lag_w015;
// Gauss-Laguerre quadrature for d=20
extern arma::vec lag_p020;
extern arma::vec lag_w020;
// Gauss-Laguerre quadrature for d=25
extern arma::vec lag_p025;
extern arma::vec lag_w025;
// Gauss-Laguerre quadrature for d=30
extern arma::vec lag_p030;
extern arma::vec lag_w030;
// Gauss-Laguerre quadrature for d=35
extern arma::vec lag_p035;
extern arma::vec lag_w035;
// Gauss-Laguerre quadrature for d=40
extern arma::vec lag_p040;
extern arma::vec lag_w040;
// Gauss-Laguerre quadrature for d=45
extern arma::vec lag_p045;
extern arma::vec lag_w045;
// Gauss-Laguerre quadrature for d=50
extern arma::vec lag_p050;
extern arma::vec lag_w050;
// Gauss-Laguerre quadrature for d=55
extern arma::vec lag_p055;
extern arma::vec lag_w055;
// Gauss-Laguerre quadrature for d=60
extern arma::vec lag_p060;
extern arma::vec lag_w060;
// Gauss-Laguerre quadrature for d=65
extern arma::vec lag_p065;
extern arma::vec lag_w065;
// Gauss-Laguerre quadrature for d=70
extern arma::vec lag_p070;
extern arma::vec lag_w070;
// Gauss-Laguerre quadrature for d=75
extern arma::vec lag_p075;
extern arma::vec lag_w075;
// Gauss-Laguerre quadrature for d=80
extern arma::vec lag_p080;
extern arma::vec lag_w080;
// Gauss-Laguerre quadrature for d=85
extern arma::vec lag_p085;
extern arma::vec lag_w085;
// Gauss-Laguerre quadrature for d=90
extern arma::vec lag_p090;
extern arma::vec lag_w090;
// Gauss-Laguerre quadrature for d=95
extern arma::vec lag_p095;
extern arma::vec lag_w095;
// Gauss-Laguerre quadrature for d=100
extern arma::vec lag_p100;
extern arma::vec lag_w100;
// Gauss-Laguerre quadrature for d=105
extern arma::vec lag_p105;
extern arma::vec lag_w105;
// Gauss-Laguerre quadrature for d=110
extern arma::vec lag_p110;
extern arma::vec lag_w110;
// Gauss-Laguerre quadrature for d=115
extern arma::vec lag_p115;
extern arma::vec lag_w115;
// Gauss-Laguerre quadrature for d=120
extern arma::vec lag_p120;
extern arma::vec lag_w120;
// Gauss-Laguerre quadrature for d=125
extern arma::vec lag_p125;
extern arma::vec lag_w125;
// Gauss-Laguerre quadrature for d=130
extern arma::vec lag_p130;
extern arma::vec lag_w130;
// Gauss-Laguerre quadrature for d=135
extern arma::vec lag_p135;
extern arma::vec lag_w135;
// Gauss-Laguerre quadrature for d=140
extern arma::vec lag_p140;
extern arma::vec lag_w140;
// Gauss-Laguerre quadrature for d=145
extern arma::vec lag_p145;
extern arma::vec lag_w145;
// Gauss-Laguerre quadrature for d=150
extern arma::vec lag_p150;
extern arma::vec lag_w150;
// Gauss-Laguerre quadrature for d=155
extern arma::vec lag_p155;
extern arma::vec lag_w155;
// Gauss-Laguerre quadrature for d=160
extern arma::vec lag_p160;
extern arma::vec lag_w160;
// Gauss-Laguerre quadrature for d=165
extern arma::vec lag_p165;
extern arma::vec lag_w165;
// Gauss-Laguerre quadrature for d=170
extern arma::vec lag_p170;
extern arma::vec lag_w170;
// Gauss-Laguerre quadrature for d=175
extern arma::vec lag_p175;
extern arma::vec lag_w175;
// Gauss-Laguerre quadrature for d=180
extern arma::vec lag_p180;
extern arma::vec lag_w180;
// Gauss-Laguerre quadrature for d=185
extern arma::vec lag_p185;
extern arma::vec lag_w185;
// Gauss-Laguerre quadrature for d=190
extern arma::vec lag_p190;
extern arma::vec lag_w190;
// Gauss-Laguerre quadrature for d=195
extern arma::vec lag_p195;
extern arma::vec lag_w195;
// Gauss-Laguerre quadrature for d=200
extern arma::vec lag_p200;
extern arma::vec lag_w200;
// Gauss-Laguerre quadrature for d=205
extern arma::vec lag_p205;
extern arma::vec lag_w205;
// Gauss-Laguerre quadrature for d=210
extern arma::vec lag_p210;
extern arma::vec lag_w210;
// Gauss-Laguerre quadrature for d=215
extern arma::vec lag_p215;
extern arma::vec lag_w215;
// Gauss-Laguerre quadrature for d=220
extern arma::vec lag_p220;
extern arma::vec lag_w220;
// Gauss-Laguerre quadrature for d=225
extern arma::vec lag_p225;
extern arma::vec lag_w225;
// Gauss-Laguerre quadrature for d=230
extern arma::vec lag_p230;
extern arma::vec lag_w230;
// Gauss-Laguerre quadrature for d=235
extern arma::vec lag_p235;
extern arma::vec lag_w235;
// Gauss-Laguerre quadrature for d=240
extern arma::vec lag_p240;
extern arma::vec lag_w240;
// Gauss-Laguerre quadrature for d=245
extern arma::vec lag_p245;
extern arma::vec lag_w245;
// Gauss-Laguerre quadrature for d=250
extern arma::vec lag_p250;
extern arma::vec lag_w250;
// Gauss-Laguerre quadrature for d=255
extern arma::vec lag_p255;
extern arma::vec lag_w255;
// Gauss-Laguerre quadrature for d=260
extern arma::vec lag_p260;
extern arma::vec lag_w260;
// Gauss-Laguerre quadrature for d=265
extern arma::vec lag_p265;
extern arma::vec lag_w265;
// Gauss-Laguerre quadrature for d=270
extern arma::vec lag_p270;
extern arma::vec lag_w270;
// Gauss-Laguerre quadrature for d=275
extern arma::vec lag_p275;
extern arma::vec lag_w275;
// Gauss-Laguerre quadrature for d=280
extern arma::vec lag_p280;
extern arma::vec lag_w280;
// Gauss-Laguerre quadrature for d=285
extern arma::vec lag_p285;
extern arma::vec lag_w285;
// Gauss-Laguerre quadrature for d=290
extern arma::vec lag_p290;
extern arma::vec lag_w290;
// Gauss-Laguerre quadrature for d=295
extern arma::vec lag_p295;
extern arma::vec lag_w295;
// Gauss-Laguerre quadrature for d=300
extern arma::vec lag_p300;
extern arma::vec lag_w300;
#endif
