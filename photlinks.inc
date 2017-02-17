  SELECT CASE (jl)
    CASE(  2)
      j(   1) = seval(n,theta,tmp,tmp2,b,c,d) ! O3 -> O2 + O(1D)
    CASE(  3)
      j(   2) = seval(n,theta,tmp,tmp2,b,c,d) ! O3 -> O2 + O(3P)
    CASE(  5)
      j(   3) = seval(n,theta,tmp,tmp2,b,c,d) ! H2O2 -> 2 OH
    CASE(  6)
      j(   4) = seval(n,theta,tmp,tmp2,b,c,d) ! NO2 -> NO + O(3P)
    CASE(  7)
      j(   5) = seval(n,theta,tmp,tmp2,b,c,d) ! NO3 -> NO + O2
    CASE(  8)
      j(   6) = seval(n,theta,tmp,tmp2,b,c,d) ! NO3 -> NO2 + O(3P)
    CASE( 12)
      j(   7) = seval(n,theta,tmp,tmp2,b,c,d) ! HNO2 -> OH + NO
    CASE( 13)
      j(   8) = seval(n,theta,tmp,tmp2,b,c,d) ! HNO3 -> OH + NO2
    CASE( 18)
      j(   9) = seval(n,theta,tmp,tmp2,b,c,d) ! HNO4 -> HO2 + NO2
    CASE( 19)
      j(  10) = seval(n,theta,tmp,tmp2,b,c,d) ! HNO4 -> OH + NO3
    CASE( 22)
      j(1001) = seval(n,theta,tmp,tmp2,b,c,d) ! CH2O -> H + HCO
    CASE( 23)
      j(1002) = seval(n,theta,tmp,tmp2,b,c,d) ! CH2O -> H2 + CO
    CASE( 24)
      j(1003) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3CHO -> CH3 + HCO
    CASE( 26)
      j(1004) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3CHO -> CH3CO + H
    CASE( 27)
      j(1005) = seval(n,theta,tmp,tmp2,b,c,d) ! C2H5CHO -> C2H5 + HCO
    CASE(131)
      j(1006) = seval(n,theta,tmp,tmp2,b,c,d) ! n-C3H7CHO -> n-C3H7 + CHO
    CASE(132)
      j(1007) = seval(n,theta,tmp,tmp2,b,c,d) ! n-C3H7CHO -> C2H4 + CH2CHOH
    CASE(140)
      j(1008) = seval(n,theta,tmp,tmp2,b,c,d) ! n-C4H9CHO -> C4H9 +  CHO
    CASE(141)
      j(1009) = seval(n,theta,tmp,tmp2,b,c,d) ! n-C4H9CHO -> CH3CH=CH2 +  CH2=CHOH
    CASE(142)
      j(1010) = seval(n,theta,tmp,tmp2,b,c,d) ! n-C4H9CHO -> 2-methylcyclobutanol
    CASE(159)
      j(1011) = seval(n,theta,tmp,tmp2,b,c,d) ! n-C5H11CHO -> C5H11 + CHO
    CASE(160)
      j(1012) = seval(n,theta,tmp,tmp2,b,c,d) ! n-C5H11CHO -> C2H5CH=CH2 + CH2=CHOH
    CASE(161)
      j(1013) = seval(n,theta,tmp,tmp2,b,c,d) ! n-C5H11CHO -> 2-ethylcyclobutanol
    CASE(167)
      j(1014) = seval(n,theta,tmp,tmp2,b,c,d) ! n-C6H13CHO -> C6H13 + CHO
    CASE(168)
      j(1015) = seval(n,theta,tmp,tmp2,b,c,d) ! n-C6H13CHO -> C3H7CH=CH2 + CH2=CHOH
    CASE(169)
      j(1016) = seval(n,theta,tmp,tmp2,b,c,d) ! n-C6H13CHO -> 2-propylcyclobutanol
    CASE(219)
      j(1017) = seval(n,theta,tmp,tmp2,b,c,d) ! nALD -> products
    CASE(137)
      j(1018) = seval(n,theta,tmp,tmp2,b,c,d) ! i-C3H7CHO -> i-C3H7 + CHO
    CASE(148)
      j(1019) = seval(n,theta,tmp,tmp2,b,c,d) ! i-C4H9CHO -> C4H9 + CHO
    CASE(149)
      j(1020) = seval(n,theta,tmp,tmp2,b,c,d) ! i-C4H9CHO -> CH3CH=CH2 + CH2=CHOH
    CASE(150)
      j(1021) = seval(n,theta,tmp,tmp2,b,c,d) ! sec-C4H9CHO -> C4H9 + CHO
    CASE(151)
      j(1022) = seval(n,theta,tmp,tmp2,b,c,d) ! sec-C4H9CHO  -> CH3CH=CHOH + CH2=CH2
    CASE(152)
      j(1023) = seval(n,theta,tmp,tmp2,b,c,d) ! t-C4H9CHO -> C4H9 + CHO
    CASE(183)
      j(1024) = seval(n,theta,tmp,tmp2,b,c,d) ! C4H9CH(C2H5)CHO -> C7H15 + CHO
    CASE(184)
      j(1025) = seval(n,theta,tmp,tmp2,b,c,d) ! C4H9CH(C2H5)CHO -> C7H16 + CO
    CASE(222)
      j(1026) = seval(n,theta,tmp,tmp2,b,c,d) ! MeALD -> NI products
    CASE(223)
      j(1027) = seval(n,theta,tmp,tmp2,b,c,d) ! MeALD -> NII products
    CASE(185)
      j(1028) = seval(n,theta,tmp,tmp2,b,c,d) ! AlkALD -> NI products
    CASE(186)
      j(1029) = seval(n,theta,tmp,tmp2,b,c,d) ! AlkALD -> NII products
    CASE(153)
      j(1030) = seval(n,theta,tmp,tmp2,b,c,d) ! tALD -> NI products
    CASE(154)
      j(1031) = seval(n,theta,tmp,tmp2,b,c,d) ! tALD -> NII products
    CASE( 49)
      j(1032) = seval(n,theta,tmp,tmp2,b,c,d) ! CH2=CHCHO -> CH2=CH + CHO
    CASE( 50)
      j(1033) = seval(n,theta,tmp,tmp2,b,c,d) ! CH2=CHCHO -> CH2=CH2 + CO
    CASE( 52)
      j(1034) = seval(n,theta,tmp,tmp2,b,c,d) ! CH2=CHCHO -> CH2=CHCO + H
    CASE(177)
      j(1035) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3CH=CHCHO -> CH3CH=CH + CHO
    CASE(178)
      j(1036) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3CH=CHCHO -> CH3CH=CH2 + CO
    CASE(179)
      j(1037) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3CH=CHCHO -> CH3CH=CHCO + H
    CASE(180)
      j(1038) = seval(n,theta,tmp,tmp2,b,c,d) ! 2-hexenal -> 1-pentenyl radical + CHO
    CASE(181)
      j(1039) = seval(n,theta,tmp,tmp2,b,c,d) ! 2-hexenal -> 1-pentene + CO
    CASE(182)
      j(1040) = seval(n,theta,tmp,tmp2,b,c,d) ! 2-hexenal -> C3H7CH=CHCO + H
    CASE(209)
      j(1041) = seval(n,theta,tmp,tmp2,b,c,d) ! hexadienal -> 1-pentenyl radical + CHO
    CASE(210)
      j(1042) = seval(n,theta,tmp,tmp2,b,c,d) ! hexadienal -> 1,3-pentadiene + CO
    CASE(211)
      j(1043) = seval(n,theta,tmp,tmp2,b,c,d) ! hexadienal -> CH3CH=CHCH=CHCO + H
    CASE( 53)
      j(1044) = seval(n,theta,tmp,tmp2,b,c,d) ! CH2=C(CH3)CHO -> CH2=CCH3 + CHO
    CASE( 54)
      j(1045) = seval(n,theta,tmp,tmp2,b,c,d) ! CH2=C(CH3)CHO -> CH3CH=CH2 + CO
    CASE( 56)
      j(1046) = seval(n,theta,tmp,tmp2,b,c,d) ! CH2=C(CH3)CHO -> CH2=C(CH3)CO + H
    CASE(200)
      j(1047) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3C(CH3)=CHCHO -> (CH3)2C=CH + CHO
    CASE(201)
      j(1048) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3C(CH3)=CHCHO -> (CH3)2C=CH2 + CO
    CASE(202)
      j(1049) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3C(CH3)=CHCHO -> (CH3)2C=CHCO + H
    CASE(191)
      j(1050) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3CH=C(CH3)CHO -> CH3CH=CCH3 + CHO
    CASE(192)
      j(1051) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3CH=C(CH3)CHO -> CH3CH=CHCH3 + CO
    CASE(193)
      j(1052) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3CH=C(CH3)CHO -> CH3CH=C(CH3)CO + H
    CASE(231)
      j(1053) = seval(n,theta,tmp,tmp2,b,c,d) ! luALD -> NI products
    CASE(232)
      j(1054) = seval(n,theta,tmp,tmp2,b,c,d) ! luALD -> alkene + CO
    CASE(233)
      j(1055) = seval(n,theta,tmp,tmp2,b,c,d) ! luALD -> acyl + H
    CASE( 61)
      j(1056) = seval(n,theta,tmp,tmp2,b,c,d) ! HOCH2CHO -> CH2OH + HCO
    CASE( 62)
      j(1057) = seval(n,theta,tmp,tmp2,b,c,d) ! HOCH2CHO -> CH3OH + CO
    CASE( 63)
      j(1058) = seval(n,theta,tmp,tmp2,b,c,d) ! HOCH2CHO -> CH2CHO + OH
    CASE(175)
      j(1059) = seval(n,theta,tmp,tmp2,b,c,d) ! Glycidaldehyde -> oxyranyl radical + CHO
    CASE(176)
      j(1060) = seval(n,theta,tmp,tmp2,b,c,d) ! Glycidaldehyde -> oxyrane + CO
    CASE( 28)
      j(1061) = seval(n,theta,tmp,tmp2,b,c,d) ! ALD3OHqy -> R(OH) + HCO
    CASE(133)
      j(1062) = seval(n,theta,tmp,tmp2,b,c,d) ! ALD4OHqy -> NI products
    CASE(134)
      j(1063) = seval(n,theta,tmp,tmp2,b,c,d) ! ALD4OHqy -> NII products
    CASE(143)
      j(1064) = seval(n,theta,tmp,tmp2,b,c,d) ! C5nALDOHqy -> NI products
    CASE(144)
      j(1065) = seval(n,theta,tmp,tmp2,b,c,d) ! C5nALDOHqy -> NII products
    CASE(162)
      j(1066) = seval(n,theta,tmp,tmp2,b,c,d) ! C6nALDOHqy -> NI products
    CASE(163)
      j(1067) = seval(n,theta,tmp,tmp2,b,c,d) ! C6nALDOHqy -> NII products
    CASE(170)
      j(1068) = seval(n,theta,tmp,tmp2,b,c,d) ! C7nALDOHqy -> NI products
    CASE(171)
      j(1069) = seval(n,theta,tmp,tmp2,b,c,d) ! C7nALDOHqy -> NII products
    CASE(174)
      j(1070) = seval(n,theta,tmp,tmp2,b,c,d) ! C7nALDOHoh -> cycl. product
    CASE(224)
      j(1071) = seval(n,theta,tmp,tmp2,b,c,d) ! MeALDOHqy -> NI products
    CASE(225)
      j(1072) = seval(n,theta,tmp,tmp2,b,c,d) ! MeALDOHqy -> NII products
    CASE(187)
      j(1073) = seval(n,theta,tmp,tmp2,b,c,d) ! AlkALDOHqy -> NI products
    CASE(188)
      j(1074) = seval(n,theta,tmp,tmp2,b,c,d) ! AlkALDOHqy -> NII products
    CASE(155)
      j(1075) = seval(n,theta,tmp,tmp2,b,c,d) ! tALDOHqy -> NI products
    CASE(156)
      j(1076) = seval(n,theta,tmp,tmp2,b,c,d) ! tALDOHqy -> NII products
    CASE(234)
      j(1077) = seval(n,theta,tmp,tmp2,b,c,d) ! luALDOHqy -> NI products
    CASE(235)
      j(1078) = seval(n,theta,tmp,tmp2,b,c,d) ! luALDOHqy -> alkene + CO
    CASE(236)
      j(1079) = seval(n,theta,tmp,tmp2,b,c,d) ! luALDOHqy -> acyl + H
    CASE(212)
      j(1080) = seval(n,theta,tmp,tmp2,b,c,d) ! uuALDOHqy -> NI products
    CASE(213)
      j(1081) = seval(n,theta,tmp,tmp2,b,c,d) ! uuALDOHqy -> diene + CO
    CASE(214)
      j(1082) = seval(n,theta,tmp,tmp2,b,c,d) ! uuALDOHqy -> acyl + H
    CASE(194)
      j(1083) = seval(n,theta,tmp,tmp2,b,c,d) ! aMeC4uALDOHqy -> NI products
    CASE(195)
      j(1084) = seval(n,theta,tmp,tmp2,b,c,d) ! aMeC4uALDOHqy -> alkene + CO
    CASE(196)
      j(1085) = seval(n,theta,tmp,tmp2,b,c,d) ! aMeC4uALDOHqy -> acyl + H
    CASE(203)
      j(1086) = seval(n,theta,tmp,tmp2,b,c,d) ! bMeC4uALDOHqy -> NI products
    CASE(204)
      j(1087) = seval(n,theta,tmp,tmp2,b,c,d) ! bMeC4uALDOHqy -> alkene + CO
    CASE(205)
      j(1088) = seval(n,theta,tmp,tmp2,b,c,d) ! bMeC4uALDOHqy -> acyl + H
    CASE(228)
      j(1089) = seval(n,theta,tmp,tmp2,b,c,d) ! cALD -> NI products
    CASE(229)
      j(1090) = seval(n,theta,tmp,tmp2,b,c,d) ! cALDOHqy -> NI products
    CASE( 70)
      j(3001) = seval(n,theta,tmp,tmp2,b,c,d) ! CHOCHO -> H2 + 2 CO
    CASE( 71)
      j(3002) = seval(n,theta,tmp,tmp2,b,c,d) ! CHOCHO -> CH2O + CO
    CASE( 69)
      j(3003) = seval(n,theta,tmp,tmp2,b,c,d) ! CHOCHO -> 2 HO2 + 2 CO
    CASE( 72)
      j(3004) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3COCHO -> CH3CO + HCO
    CASE( 73)
      j(3005) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3COCOCH3 -> Products
    CASE(307)
      j(3006) = seval(n,theta,tmp,tmp2,b,c,d) ! CHOCH=CHCHO -> 3H-furan-2-one
    CASE(308)
      j(3007) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3COCH=CHCHO -> 5Me-3H-2-furanone
    CASE(309)
      j(3008) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3COCH=CHCHO -> CH3 + CHOCH=CHCO
    CASE(310)
      j(3009) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3COCH=CHCHO -> CH3COCH=CH2 + CO
    CASE(313)
      j(3010) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3COCH=CHCOCH3 -> CH3CO + CH=CHCOCH3
    CASE(312)
      j(3011) = seval(n,theta,tmp,tmp2,b,c,d) ! CHOCH=CHCH=CHCHO -> diformyl cyclobutene
    CASE(138)
      j(3012) = seval(n,theta,tmp,tmp2,b,c,d) ! pinonaldehyde -> R + CO + HO2
    CASE(139)
      j(3013) = seval(n,theta,tmp,tmp2,b,c,d) ! caronaldehyde -> R + CO + HO2
    CASE( 32)
      j(4001) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3ONO2 -> CH3O + NO2
    CASE( 36)
      j(4002) = seval(n,theta,tmp,tmp2,b,c,d) ! C2H5ONO2 -> C2H5O + NO2
    CASE( 37)
      j(4003) = seval(n,theta,tmp,tmp2,b,c,d) ! n-C3H7ONO2 -> C3H7O + NO2
    CASE( 38)
      j(4004) = seval(n,theta,tmp,tmp2,b,c,d) ! 1-C4H9ONO2 -> 1-C4H9O + NO2
    CASE(314)
      j(4005) = seval(n,theta,tmp,tmp2,b,c,d) ! n-C5H11ONO2 -> n-C5H11O + NO2
    CASE( 40)
      j(4006) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3CHONO2CH3 -> CH3CHOCH3 + NO2
    CASE( 39)
      j(4007) = seval(n,theta,tmp,tmp2,b,c,d) ! 2-C4H9ONO2 -> 2-C4H9O + NO2
    CASE(315)
      j(4008) = seval(n,theta,tmp,tmp2,b,c,d) ! 2-C5H11ONO2 -> 2-C5H11O + NO2
    CASE(316)
      j(4009) = seval(n,theta,tmp,tmp2,b,c,d) ! 3-C5H11ONO2 -> 3-C5H11O + NO2
    CASE(318)
      j(4010) = seval(n,theta,tmp,tmp2,b,c,d) ! i-C4H9ONO2 -> i-C4H9O + NO2
    CASE( 43)
      j(4011) = seval(n,theta,tmp,tmp2,b,c,d) ! C(CH3)3(ONO2) -> C(CH3)3(O.) + NO2
    CASE(319)
      j(4012) = seval(n,theta,tmp,tmp2,b,c,d) ! i-C5H11ONO2 -> i-C5H11O + NO2
    CASE(317)
      j(4013) = seval(n,theta,tmp,tmp2,b,c,d) ! c-C5H11ONO2 -> c-C5H11O + NO2
    CASE( 41)
      j(4014) = seval(n,theta,tmp,tmp2,b,c,d) ! CH2(OH)CH2(ONO2) -> CH2(OH)CH2(O.) + NO2
    CASE(320)
      j(4015) = seval(n,theta,tmp,tmp2,b,c,d) ! C1(OH)NO3 -> C1(OH)O + NO2
    CASE(321)
      j(4016) = seval(n,theta,tmp,tmp2,b,c,d) ! R(OH)NO3 -> R(OH)O + NO2
    CASE(322)
      j(4017) = seval(n,theta,tmp,tmp2,b,c,d) ! iR(OH)NO3 -> iR(OH)O + NO2
    CASE(323)
      j(4018) = seval(n,theta,tmp,tmp2,b,c,d) ! tR(OH)NO3 -> tR(OH)O + NO2
    CASE( 35)
      j(4019) = seval(n,theta,tmp,tmp2,b,c,d) ! ROONO2 -> products
    CASE( 45)
      j(4020) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3CO(OONO2) -> CH3CO(OO) + NO2
    CASE( 46)
      j(4021) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3CO(OONO2) -> CH3CO(O) + NO3
    CASE( 47)
      j(4022) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3CH2CO(OONO2) -> CH3CH2CO(OO) + NO2
    CASE( 48)
      j(4023) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3CH2CO(OONO2) -> CH3CH2CO(O) + NO3
    CASE(324)
      j(4024) = seval(n,theta,tmp,tmp2,b,c,d) ! PAN -> RCO(OO) + NO2
    CASE(325)
      j(4025) = seval(n,theta,tmp,tmp2,b,c,d) ! PAN -> RCO(O) + NO3
    CASE(326)
      j(5001) = seval(n,theta,tmp,tmp2,b,c,d)*2.000 ! CH3CH(NO3)CH2NO3 -> CH3CH(NO3)CH2O + NO2
    CASE(328)
      j(5002) = seval(n,theta,tmp,tmp2,b,c,d)*2.000 ! C2H5CH(NO3)CH2NO3 -> C2H5CH(NO3)CH2O + NO2
    CASE(330)
      j(5003) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3CH(NO3)CH(NO3)CH3 -> RO. + NO2
    CASE(335)
      j(5004) = seval(n,theta,tmp,tmp2,b,c,d)*2.000 ! C6H9-1-CH3-1,2-NO3 -> R1O. + NO2
    CASE(331)
      j(5005) = seval(n,theta,tmp,tmp2,b,c,d) ! CH2(NO3)CH=CHCH2NO3 -> RO. + NO2
    CASE(332)
      j(5006) = seval(n,theta,tmp,tmp2,b,c,d)*2.000 ! CH2=CHCH(NO3)CH2NO3 -> C2H3CH(NO3)CH2O + NO2
    CASE(334)
      j(5007) = seval(n,theta,tmp,tmp2,b,c,d) ! uDINIT -> RO. + NO2
    CASE( 30)
      j(6001) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3OOH -> CH3O + OH
    CASE(337)
      j(6002) = seval(n,theta,tmp,tmp2,b,c,d) ! (CH3)3COOH -> (CH3)3CO + OH
    CASE( 31)
      j(6003) = seval(n,theta,tmp,tmp2,b,c,d) ! HOCH2OOH -> HOCH2O. + OH
    CASE( 75)
      j(6004) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3CO(OOH) -> Products
    CASE(343)
      j(7001) = seval(n,theta,tmp,tmp2,b,c,d) ! CH2OO -> HCHO + O(3P)
    CASE(344)
      j(7002) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3CHOO -> CH3CHO + O(3P)
    CASE(345)
      j(7004) = seval(n,theta,tmp,tmp2,b,c,d) ! C2H5CHOO -> C2H5CHO + O(3P)
    CASE(346)
      j(7005) = seval(n,theta,tmp,tmp2,b,c,d) ! (CH3)2COO -> CH3COCH3 + O(3P)
    CASE( 57)
      j(  20) = seval(n,theta,tmp,tmp2,b,c,d) ! HPALD -> Products
  END SELECT
