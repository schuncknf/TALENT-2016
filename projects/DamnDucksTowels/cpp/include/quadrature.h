#ifndef QUADRATURES_H
#define QUADRATURES_H

#define GET_LAG_ROOTS(d, xipp, wipp) switch(d) { case 5: xipp = &lag_p005; wipp = &lag_w005; break; case 10: xipp = &lag_p010; wipp = &lag_w010; break; case 15: xipp = &lag_p015; wipp = &lag_w015; break; case 20: xipp = &lag_p020; wipp = &lag_w020; break; case 25: xipp = &lag_p025; wipp = &lag_w025; break; case 30: xipp = &lag_p030; wipp = &lag_w030; break; case 35: xipp = &lag_p035; wipp = &lag_w035; break; case 40: xipp = &lag_p040; wipp = &lag_w040; break; case 45: xipp = &lag_p045; wipp = &lag_w045; break; case 50: xipp = &lag_p050; wipp = &lag_w050; break; case 55: xipp = &lag_p055; wipp = &lag_w055; break; case 60: xipp = &lag_p060; wipp = &lag_w060; break; case 65: xipp = &lag_p065; wipp = &lag_w065; break; case 70: xipp = &lag_p070; wipp = &lag_w070; break; case 75: xipp = &lag_p075; wipp = &lag_w075; break; case 80: xipp = &lag_p080; wipp = &lag_w080; break; case 85: xipp = &lag_p085; wipp = &lag_w085; break; case 90: xipp = &lag_p090; wipp = &lag_w090; break; case 95: xipp = &lag_p095; wipp = &lag_w095; break; case 100: xipp = &lag_p100; wipp = &lag_w100; break; case 105: xipp = &lag_p105; wipp = &lag_w105; break; case 110: xipp = &lag_p110; wipp = &lag_w110; break; case 115: xipp = &lag_p115; wipp = &lag_w115; break; case 120: xipp = &lag_p120; wipp = &lag_w120; break; case 125: xipp = &lag_p125; wipp = &lag_w125; break; case 130: xipp = &lag_p130; wipp = &lag_w130; break; case 135: xipp = &lag_p135; wipp = &lag_w135; break; case 140: xipp = &lag_p140; wipp = &lag_w140; break; case 145: xipp = &lag_p145; wipp = &lag_w145; break; case 150: xipp = &lag_p150; wipp = &lag_w150; break; case 155: xipp = &lag_p155; wipp = &lag_w155; break; case 160: xipp = &lag_p160; wipp = &lag_w160; break; case 165: xipp = &lag_p165; wipp = &lag_w165; break; case 170: xipp = &lag_p170; wipp = &lag_w170; break; case 175: xipp = &lag_p175; wipp = &lag_w175; break; case 180: xipp = &lag_p180; wipp = &lag_w180; break; case 185: xipp = &lag_p185; wipp = &lag_w185; break; case 190: xipp = &lag_p190; wipp = &lag_w190; break; case 195: xipp = &lag_p195; wipp = &lag_w195; break; case 200: xipp = &lag_p200; wipp = &lag_w200; break; default: xipp = NULL; wipp = NULL; }
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
#endif
