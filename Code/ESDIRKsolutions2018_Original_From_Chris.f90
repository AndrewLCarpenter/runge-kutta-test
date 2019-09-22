c
c ===============================================
c THIRD-ORDER
c ===============================================
c
c ESDIRK3(2)4L[2]SA 
c
      b(1) = 1471266399579.d0/7840856788654.d0
      b(2) = -4482444167858.d0/7529755066697.d0
      b(3) = 11266239266428.d0/11593286722821.d0
      b(4) = 1767732205903.d0/4055673282236.d0
      c(1) = 0.d0/1.d0
      c(2) = 1767732205903.d0/2027836641118.d0
      c(3) = 3.d0/5.d0
      c(4) = 1.d0/1.d0
c
      bh(1) = 926040629867.d0/8503851176844.d0
      bh(2) = -19534562426408.d0/21341649249991.d0
      bh(3) = 17036650473653.d0/13401246206802.d0
      bh(4) = 4543788980243.d0/8490594148910.d0
c
      a(1,1) = 0.d0/1.d0
      a(2,1) = 1767732205903.d0/4055673282236.d0
      a(2,2) = 1767732205903.d0/4055673282236.d0
      a(3,1) = 2746238789719.d0/10658868560708.d0
      a(3,2) = -640167445237.d0/6845629431997.d0
      a(3,3) = 1767732205903.d0/4055673282236.d0
      a(4,1) = 1471266399579.d0/7840856788654.d0
      a(4,2) = -4482444167858.d0/7529755066697.d0
      a(4,3) = 11266239266428.d0/11593286722821.d0
      a(4,4) = 1767732205903.d0/4055673282236.d0
c
c ESDIRK3(2)5L[2]SA (page 85, TM-2016-219173)
c
      b(1) = 1150014452202.d0/6551090391107.d0
      b(2) = 1150014452202.d0/6551090391107.d0
      b(3) = -1642152774673.d0/4734363421565.d0
      b(4) = 5827.d0/7560.d0
      b(5) = 9.d0/40.d0
      c(1) = 0.d0/1.d0
      c(2) = 9.d0/20.d0
      c(3) = 5627749853384.d0/7325910085487.d0
      c(4) = 3.d0/5.d0
      c(5) = 1.d0/1.d0
      bh(1) = 4555948517383.d0/24713416420891.d0
      bh(2) = 4555948517383.d0/24713416420891.d0
      bh(3) = -7107561914881.d0/25547637784726.d0
      bh(4) = 30698249.d0/44052120.d0
      bh(5) = 49563.d0/233080.d0
      a(1,1) = 0.d0/1.d0
      a(2,1) = 9.d0/40.d0
      a(2,2) = 9.d0/40.d0
      a(3,1) = 1748874742213.d0/6439178996590.d0
      a(3,2) = 1748874742213.d0/6439178996590.d0
      a(3,3) = 9.d0/40.d0
      a(4,1) = 2511125995703.d0/11223226150663.d0
      a(4,2) = 2511125995703.d0/11223226150663.d0
      a(4,3) = -683874557734.d0/9434395612819.d0
      a(4,4) = 9.d0/40.d0
      a(5,1) = 1150014452202.d0/6551090391107.d0
      a(5,2) = 1150014452202.d0/6551090391107.d0
      a(5,3) = -1642152774673.d0/4734363421565.d0
      a(5,4) = 5827.d0/7560.d0
      a(5,5) = 9.d0/40.d0
c
c ===============================================
c FOURTH-ORDER
c ===============================================
c
c SDIRK4 - Hairer & Wanner(1996)
c
      b(1) = 25.d0/24.d0
      b(2) = -49.d0/48.d0
      b(3) = 125.d0/16.d0
      b(4) = -85.d0/12.d0
      b(5) = 1.d0/4.d0
      c(1) = 1.d0/4.d0
      c(2) = 3.d0/4.d0
      c(3) = 11.d0/20.d0
      c(4) = 1.d0/2.d0
      c(5) = 1.d0/1.d0
      bh(1) = 59.d0/48.d0
      bh(2) = -17.d0/96.d0
      bh(3) = 225.d0/32.d0
      bh(4) = -85.d0/12.d0
      bh(5) = 0.d0/1.d0
      a(1,1) = 1.d0/4.d0
      a(2,1) = 1.d0/2.d0
      a(2,2) = 1.d0/4.d0
      a(3,1) = 17.d0/50.d0
      a(3,2) = -1.d0/25.d0
      a(3,3) = 1.d0/4.d0
      a(4,1) = 371.d0/1360.d0
      a(4,2) = -137.d0/2720.d0
      a(4,3) = 15.d0/544.d0
      a(4,4) = 1.d0/4.d0
      a(5,1) = 25.d0/24.d0
      a(5,2) = -49.d0/48.d0
      a(5,3) = 125.d0/16.d0
      a(5,4) = -85.d0/12.d0
      a(5,5) = 1.d0/4.d0
c
c ESDIRK4(3)6L[2]SA_ARK
c
c     b(1) = 82889.d0/524892.d0
c     b(2) = 0.d0/1.d0
c     b(3) = 15625.d0/83664.d0
c     b(4) = 69875.d0/102672.d0
c     b(5) = -2260.d0/8211.d0
c     b(6) = 1.d0/4.d0
c     c(1) = 0.d0/1.d0
c     c(2) = 1.d0/2.d0
c     c(3) = 83.d0/250.d0
c     c(4) = 31.d0/50.d0
c     c(5) = 17.d0/20.d0
c     c(6) = 1.d0/1.d0
c     bh(1) = 4586570599.d0/29645900160.d0
c     bh(2) = 0.d0/1.d0
c     bh(3) = 178811875.d0/945068544.d0
c     bh(4) = 814220225.d0/1159782912.d0
c     bh(5) = -3700637.d0/11593932.d0
c     bh(6) = 61727.d0/225920.d0
c     a(1,1) = 0.d0/1.d0
c     a(2,1) = 1.d0/4.d0
c     a(2,2) = 1.d0/4.d0
c     a(3,1) = 8611.d0/62500.d0
c     a(3,2) = -1743.d0/31250.d0
c     a(3,3) = 1.d0/4.d0
c     a(4,1) = 5012029.d0/34652500.d0
c     a(4,2) = -654441.d0/2922500.d0
c     a(4,3) = 174375.d0/388108.d0
c     a(4,4) = 1.d0/4.d0
c     a(5,1) = 15267082809.d0/155376265600.d0
c     a(5,2) = -71443401.d0/120774400.d0
c     a(5,3) = 730878875.d0/902184768.d0
c     a(5,4) = 2285395.d0/8070912.d0
c     a(5,5) = 1.d0/4.d0
c     a(6,1) = 82889.d0/524892.d0
c     a(6,2) = 0.d0/1.d0
c     a(6,3) = 15625.d0/83664.d0
c     a(6,4) = 69875.d0/102672.d0
c     a(6,5) = -2260.d0/8211.d0
c     a(6,6) = 1.d0/4.d0
c
c ESDIRK4(3)6L[2]SA_1 (page 90, TM-2016-219173)
c
      b(1) = -140404485182.d0/9007427031765.d0
      b(2) = -140404485182.d0/9007427031765.d0
      b(3) = 1402990318619.d0/3619147572429.d0
      b(4) = 1394850835967.d0/2779846451479.d0
      b(5) = -990969975725.d0/9154032505244.d0
      b(6) = 1.d0/4.d0
      c(1) = 0.d0/1.d0
      c(2) = 1.d0/2.d0
      c(3) = 1513744654945.d0/10336495061766.d0
      c(4) = 5.d0/8.d0
      c(5) = 26.d0/25.d0
      c(6) = 1.d0/1.d0
      bh(1) = -140404485182.d0/9007427031765.d0
      bh(2) = -140404485182.d0/9007427031765.d0
      bh(3) = 1402990318619.d0/3619147572429.d0
      bh(4) = 1394850835967.d0/2779846451479.d0
      bh(5) = -990969975725.d0/9154032505244.d0
      bh(6) = 1.d0/4.d0
      a(1,1) = 0.d0/1.d0
      a(2,1) = 1.d0/4.d0
      a(2,2) = 1.d0/4.d0
      a(3,1) = -1356991263433.d0/26208533697614.d0
      a(3,2) = -1356991263433.d0/26208533697614.d0
      a(3,3) = 1.d0/4.d0
      a(4,1) = -162708388469.d0/2125389860943.d0
      a(4,2) = -162708388469.d0/2125389860943.d0
      a(4,3) = 1768468916152.d0/3348680272939.d0
      a(4,4) = 1.d0/4.d0
      a(5,1) = -1410042975763.d0/1938452943086.d0
      a(5,2) = -1410042975763.d0/1938452943086.d0
      a(5,3) = 20915361222208.d0/13195852609937.d0
      a(5,4) = 7240912463611.d0/10974111771892.d0
      a(5,5) = 1.d0/4.d0
      a(6,1) = -140404485182.d0/9007427031765.d0
      a(6,2) = -140404485182.d0/9007427031765.d0
      a(6,3) = 1402990318619.d0/3619147572429.d0
      a(6,4) = 1394850835967.d0/2779846451479.d0
      a(6,5) = -990969975725.d0/9154032505244.d0
      a(6,6) = 1.d0/4.d0
c
c ESDIRK4(3)6L[2]SA_2 (2018 paper)
c
      b(1) = -12917657251.d0/5222094901039.d0
      b(2) = -12917657251.d0/5222094901039.d0
      b(3) = 5602338284630.d0/15643096342197.d0
      b(4) = 9002339615474.d0/18125249312447.d0
      b(5) = -2420307481369.d0/24731958684496.d0
      b(6) = 31.d0/125.d0
      c(1) = 0.d0/1.d0
      c(2) = 62.d0/125.d0
      c(3) = 486119545908.d0/3346201505189.d0
      c(4) = 1043.d0/1706.d0
      c(5) = 1361.d0/1300.d0
      c(6) = 1.d0/1.d0
      bh(1) = -1007911106287.d0/12117826057527.d0
      bh(2) = -1007911106287.d0/12117826057527.d0
      bh(3) = 17694008993113.d0/35931961998873.d0
      bh(4) = 5816803040497.d0/11256217655929.d0
      bh(5) = -538664890905.d0/7490061179786.d0
      bh(6) = 2032560730450.d0/8872919773257.d0
      a(1,1) = 0.d0/1.d0
      a(2,1) = 31.d0/125.d0
      a(2,2) = 31.d0/125.d0
      a(3,1) = -360286518617.d0/7014585480527.d0
      a(3,2) = -360286518617.d0/7014585480527.d0
      a(3,3) = 31.d0/125.d0
      a(4,1) = -506388693497.d0/5937754990171.d0
      a(4,2) = -506388693497.d0/5937754990171.d0
      a(4,3) = 7149918333491.d0/13390931526268.d0
      a(4,4) = 31.d0/125.d0
      a(5,1) = -7628305438933.d0/11061539393788.d0
      a(5,2) = -7628305438933.d0/11061539393788.d0
      a(5,3) = 21592626537567.d0/14352247503901.d0
      a(5,4) = 11630056083252.d0/17263101053231.d0
      a(5,5) = 31.d0/125.d0
      a(6,1) = -12917657251.d0/5222094901039.d0
      a(6,2) = -12917657251.d0/5222094901039.d0
      a(6,3) = 5602338284630.d0/15643096342197.d0
      a(6,4) = 9002339615474.d0/18125249312447.d0
      a(6,5) = -2420307481369.d0/24731958684496.d0
      a(6,6) = 31.d0/125.d0
c
c ESDIRK4(3)7L[2]SA_1 (2018 paper)
c
      b(1) = -5649241495537.d0/14093099002237.d0
      b(2) = -5649241495537.d0/14093099002237.d0
      b(3) = 5718691255176.d0/6089204655961.d0
      b(4) = 2199600963556.d0/4241893152925.d0
      b(5) = 8860614275765.d0/11425531467341.d0
      b(6) = -3696041814078.d0/6641566663007.d0
      b(7) = 1.d0/8.d0
      c(1) = 0.d0/1.d0
      c(2) = 1.d0/4.d0
      c(3) = 1200237871921.d0/16391473681546.d0
      c(4) = 1.d0/2.d0
      c(5) = 395.d0/567.d0
      c(6) = 89.d0/126.d0
      c(7) = 1.d0/1.d0
      bh(1) = -1517409284625.d0/6267517876163.d0
      bh(2) = -1517409284625.d0/6267517876163.d0
      bh(3) = 8291371032348.d0/12587291883523.d0
      bh(4) = 5328310281212.d0/10646448185159.d0
      bh(5) = 5405006853541.d0/7104492075037.d0
      bh(6) = -4254786582061.d0/7445269677723.d0
      bh(7) = 19.d0/140.d0
      a(1,1) = 0.d0/1.d0
      a(2,1) = 1.d0/8.d0
      a(2,2) = 1.d0/8.d0
      a(3,1) = -39188347878.d0/1513744654945.d0
      a(3,2) = -39188347878.d0/1513744654945.d0
      a(3,3) = 1.d0/8.d0
      a(4,1) = 1748874742213.d0/5168247530883.d0
      a(4,2) = 1748874742213.d0/5168247530883.d0
      a(4,3) = -1748874742213.d0/5795261096931.d0
      a(4,4) = 1.d0/8.d0
      a(5,1) = -6429340993097.d0/17896796106705.d0
      a(5,2) = -6429340993097.d0/17896796106705.d0
      a(5,3) = 9711656375562.d0/10370074603625.d0
      a(5,4) = 1137589605079.d0/3216875020685.d0
      a(5,5) = 1.d0/8.d0
      a(6,1) = 405169606099.d0/1734380148729.d0
      a(6,2) = 405169606099.d0/1734380148729.d0
      a(6,3) = -264468840649.d0/6105657584947.d0
      a(6,4) = 118647369377.d0/6233854714037.d0
      a(6,5) = 683008737625.d0/4934655825458.d0
      a(6,6) = 1.d0/8.d0
      a(7,1) = -5649241495537.d0/14093099002237.d0
      a(7,2) = -5649241495537.d0/14093099002237.d0
      a(7,3) = 5718691255176.d0/6089204655961.d0
      a(7,4) = 2199600963556.d0/4241893152925.d0
      a(7,5) = 8860614275765.d0/11425531467341.d0
      a(7,6) = -3696041814078.d0/6641566663007.d0
      a(7,7) = 1.d0/8.d0
c
c ===============================================
c FIFTH-ORDER
c ===============================================
c
c ESDIRK5(4)7L[2]SA_1 (2018 paper)
c
      b(1) = -1319096626979.d0/17356965168099.d0
      b(2) = -1319096626979.d0/17356965168099.d0
      b(3) = 4356877330928.d0/10268933656267.d0
      b(4) = 922991294344.d0/3350617878647.d0
      b(5) = 4729382008034.d0/14755765856909.d0
      b(6) = -308199069217.d0/5897303561678.d0
      b(7) = 23.d0/125.d0
      bh(1) = -1411771963762.d0/13065992744575.d0
      bh(2) = -1411771963762.d0/13065992744575.d0
      bh(3) = 3444147136021.d0/7120014004471.d0
      bh(4) = 513148406631.d0/2174808674020.d0
      bh(5) = 2235823584917.d0/5956107268851.d0
      bh(6) = -409370279021.d0/12671388721973.d0
      bh(7) = 1600355288163.d0/10436417755451.d0
      c(1) = 0.d0/1.d0
      c(2) = 46.d0/125.d0
      c(3) = 1518047795759.d0/14084074382095.d0
      c(4) = 13.d0/25.d0
      c(5) = 5906118540659.d0/9042400211275.d0
      c(6) = 26.d0/25.d0
      c(7) = 1.d0/1.d0
      a(1,1) = 0.d0/1.d0
      a(2,1) = 23.d0/125.d0
      a(2,2) = 23.d0/125.d0
      a(3,1) = -121529886477.d0/3189120653983.d0
      a(3,2) = -121529886477.d0/3189120653983.d0
      a(3,3) = 23.d0/125.d0
      a(4,1) = 186345625210.d0/8596203768457.d0
      a(4,2) = 186345625210.d0/8596203768457.d0
      a(4,3) = 3681435451073.d0/12579882114497.d0
      a(4,4) = 23.d0/125.d0
      a(5,1) = -9898129553915.d0/11630542248213.d0
      a(5,2) = -9898129553915.d0/11630542248213.d0
      a(5,3) = 19565727496993.d0/11159348038501.d0
      a(5,4) = 2073446517052.d0/4961027473423.d0
      a(5,5) = 23.d0/125.d0
      a(6,1) = -39752543191591.d0/7894275939720.d0
      a(6,2) = -39752543191591.d0/7894275939720.d0
      a(6,3) = 52228808998390.d0/5821762529307.d0
      a(6,4) = 2756378382725.d0/8748785577174.d0
      a(6,5) = 17322065038796.d0/10556643942083.d0
      a(6,6) = 23.d0/125.d0
      a(7,1) = -1319096626979.d0/17356965168099.d0
      a(7,2) = -1319096626979.d0/17356965168099.d0
      a(7,3) = 4356877330928.d0/10268933656267.d0
      a(7,4) = 922991294344.d0/3350617878647.d0
      a(7,5) = 4729382008034.d0/14755765856909.d0
      a(7,6) = -308199069217.d0/5897303561678.d0
      a(7,7) = 23.d0/125.d0
c
c ESDIRK5(4)7L[2]SA_2 (2018 paper)
c
      b(1) = -188593204321.d0/4778616380481.d0
      b(2) = -188593204321.d0/4778616380481.d0
      b(3) = 2809310203510.d0/10304234040467.d0
      b(4) = 1021729336898.d0/2364210264653.d0
      b(5) = 870612361811.d0/2470410392208.d0
      b(6) = -1307970675534.d0/8059683598661.d0
      b(7) = 23.d0/125.d0
      bh(1) = -582099335757.d0/7214068459310.d0
      bh(2) = -582099335757.d0/7214068459310.d0
      bh(3) = 615023338567.d0/3362626566945.d0
      bh(4) = 3192122436311.d0/6174152374399.d0
      bh(5) = 6156034052041.d0/14430468657929.d0
      bh(6) = -1011318518279.d0/9693750372484.d0
      bh(7) = 1914490192573.d0/13754262428401.d0
      c(1) = 0.d0/1.d0
      c(2) = 46.d0/125.d0
      c(3) = 7121331996143.d0/11335814405378.d0
      c(4) = 49.d0/353.d0
      c(5) = 3706679970760.d0/5295570149437.d0
      c(6) = 347.d0/382.d0
      c(7) = 1.d0/1.d0
      a(1,1) = 0.d0/1.d0
      a(2,1) = 23.d0/125.d0
      a(2,2) = 23.d0/125.d0
      a(3,1) = 791020047304.d0/3561426431547.d0
      a(3,2) = 791020047304.d0/3561426431547.d0
      a(3,3) = 23.d0/125.d0
      a(4,1) = -158159076358.d0/11257294102345.d0
      a(4,2) = -158159076358.d0/11257294102345.d0
      a(4,3) = -85517644447.d0/5003708988389.d0
      a(4,4) = 23.d0/125.d0
      a(5,1) = -1653327111580.d0/4048416487981.d0
      a(5,2) = -1653327111580.d0/4048416487981.d0
      a(5,3) = 1514767744496.d0/9099671765375.d0
      a(5,4) = 14283835447591.d0/12247432691556.d0
      a(5,5) = 23.d0/125.d0
      a(6,1) = -4540011970825.d0/8418487046959.d0
      a(6,2) = -4540011970825.d0/8418487046959.d0
      a(6,3) = -1790937573418.d0/7393406387169.d0
      a(6,4) = 10819093665085.d0/7266595846747.d0
      a(6,5) = 4109463131231.d0/7386972500302.d0
      a(6,6) = 23.d0/125.d0
      a(7,1) = -188593204321.d0/4778616380481.d0
      a(7,2) = -188593204321.d0/4778616380481.d0
      a(7,3) = 2809310203510.d0/10304234040467.d0
      a(7,4) = 1021729336898.d0/2364210264653.d0
      a(7,5) = 870612361811.d0/2470410392208.d0
      a(7,6) = -1307970675534.d0/8059683598661.d0
      a(7,7) = 23.d0/125.d0
c
c ESDIRK5(4)8L[2]SA (2018 paper)
c
      b(1) = 2162042939093.d0/22873479087181.d0
      b(2) = 2162042939093.d0/22873479087181.d0
      b(3) = -4222515349147.d0/9397994281350.d0
      b(4) = 3431955516634.d0/4748630552535.d0
      b(5) = -374165068070.d0/9085231819471.d0
      b(6) = -1847934966618.d0/8254951855109.d0
      b(7) = 5186241678079.d0/7861334770480.d0
      b(8) = 1.d0/7.d0
      c(1) = 0.d0/1.d0
      c(2) = 2.d0/7.d0
      c(3) = (2.d0 + sqrt(2.d0))/7.d0
      c(4) = 150.d0/203.d0
      c(5) = 27.d0/46.d0
      c(6) = 473.d0/532.d0
      c(7) = 30.d0/83.d0
      c(8) = 1.d0/1.d0
      bh(1) = 701879993119.d0/7084679725724.d0
      bh(2) = 701879993119.d0/7084679725724.d0
      bh(3) = -8461269287478.d0/14654112271769.d0
      bh(4) = 6612459227430.d0/11388259134383.d0
      bh(5) = 2632441606103.d0/12598871370240.d0
      bh(6) = -2147694411931.d0/10286892713802.d0
      bh(7) = 4103061625716.d0/6371697724583.d0
      bh(8) = 36.d0/233.d0
      a(1,1) = 0.d0/1.d0
      a(2,1) = 1.d0/7.d0
      a(2,2) = 1.d0/7.d0
      a(3,1) = 1521428834970.d0/8822750406821.d0
      a(3,2) = 1521428834970.d0/8822750406821.d0
      a(3,3) = 1.d0/7.d0
      a(4,1) = 5338711108027.d0/29869763600956.d0
      a(4,2) = 5338711108027.d0/29869763600956.d0
      a(4,3) = 1483184435021.d0/6216373359362.d0
      a(4,4) = 1.d0/7.d0
      a(5,1) = 2264935805846.d0/12599242299355.d0
      a(5,2) = 2264935805846.d0/12599242299355.d0
      a(5,3) = 1330937762090.d0/13140498839569.d0
      a(5,4) = -287786842865.d0/17211061626069.d0
      a(5,5) = 1.d0/7.d0
      a(6,1) = 118352937080.d0/527276862197.d0
      a(6,2) = 118352937080.d0/527276862197.d0
      a(6,3) = -2960446233093.d0/7419588050389.d0
      a(6,4) = -3064256220847.d0/46575910191280.d0
      a(6,5) = 6010467311487.d0/7886573591137.d0
      a(6,6) = 1.d0/7.d0
      a(7,1) = 1134270183919.d0/9703695183946.d0
      a(7,2) = 1134270183919.d0/9703695183946.d0
      a(7,3) = 4862384331311.d0/10104465681802.d0
      a(7,4) = 1127469817207.d0/2459314315538.d0
      a(7,5) = -9518066423555.d0/11243131997224.d0
      a(7,6) = -811155580665.d0/7490894181109.d0
      a(7,7) = 1.d0/7.d0
      a(8,1) = 2162042939093.d0/22873479087181.d0
      a(8,2) = 2162042939093.d0/22873479087181.d0
      a(8,3) = -4222515349147.d0/9397994281350.d0
      a(8,4) = 3431955516634.d0/4748630552535.d0
      a(8,5) = -374165068070.d0/9085231819471.d0
      a(8,6) = -1847934966618.d0/8254951855109.d0
      a(8,7) = 5186241678079.d0/7861334770480.d0
      a(8,8) = 1.d0/7.d0
