  SELECT CASE (jl)
    CASE(  2)
      j(   1) = seval(szabin,theta,tmp,tmp2,b,c,d) ! O3 -> O2 + O(1D)
    CASE(  3)
      j(   2) = seval(szabin,theta,tmp,tmp2,b,c,d) ! O3 -> O2 + O(3P)
    CASE(  5)
      j(   3) = seval(szabin,theta,tmp,tmp2,b,c,d) ! H2O2 -> 2 OH
    CASE(  6)
      j(   4) = seval(szabin,theta,tmp,tmp2,b,c,d) ! NO2 -> NO + O(3P)
    CASE(  7)
      j(   5) = seval(szabin,theta,tmp,tmp2,b,c,d) ! NO3 -> NO + O2
    CASE(  8)
      j(   6) = seval(szabin,theta,tmp,tmp2,b,c,d) ! NO3 -> NO2 + O(3P)
    CASE( 12)
      j(   7) = seval(szabin,theta,tmp,tmp2,b,c,d) ! HNO2 -> OH + NO
    CASE( 13)
      j(   8) = seval(szabin,theta,tmp,tmp2,b,c,d) ! HNO3 -> OH + NO2
    CASE( 18)
      j(   9) = seval(szabin,theta,tmp,tmp2,b,c,d) ! HNO4 -> HO2 + NO2
    CASE( 19)
      j(  10) = seval(szabin,theta,tmp,tmp2,b,c,d) ! HNO4 -> OH + NO3
    CASE( 22)
      j(1001) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH2O -> H + HCO
    CASE( 23)
      j(1002) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH2O -> H2 + CO
    CASE( 24)
      j(1003) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3CHO -> CH3 + HCO
    CASE( 26)
      j(1004) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3CHO -> CH3CO + H
    CASE( 28)
      j(1005) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C2H5CHO -> C2H5 + HCO
    CASE(135)
      j(1006) = seval(szabin,theta,tmp,tmp2,b,c,d) ! n-C3H7CHO -> n-C3H7 + CHO
    CASE(136)
      j(1007) = seval(szabin,theta,tmp,tmp2,b,c,d) ! n-C3H7CHO -> C2H4 + CH2CHOH
    CASE(144)
      j(1008) = seval(szabin,theta,tmp,tmp2,b,c,d) ! n-C4H9CHO -> C4H9 +  CHO
    CASE(145)
      j(1009) = seval(szabin,theta,tmp,tmp2,b,c,d) ! n-C4H9CHO -> CH3CH=CH2 +  CH2=CHOH
    CASE(146)
      j(1010) = seval(szabin,theta,tmp,tmp2,b,c,d) ! n-C4H9CHO -> 2-methylcyclobutanol
    CASE(160)
      j(1011) = seval(szabin,theta,tmp,tmp2,b,c,d) ! n-C5H11CHO -> C5H11 + CHO
    CASE(161)
      j(1012) = seval(szabin,theta,tmp,tmp2,b,c,d) ! n-C5H11CHO -> C2H5CH=CH2 + CH2=CHOH
    CASE(162)
      j(1013) = seval(szabin,theta,tmp,tmp2,b,c,d) ! n-C5H11CHO -> 2-ethylcyclobutanol
    CASE(169)
      j(1014) = seval(szabin,theta,tmp,tmp2,b,c,d) ! n-C6H13CHO -> C6H13 + CHO
    CASE(170)
      j(1015) = seval(szabin,theta,tmp,tmp2,b,c,d) ! n-C6H13CHO -> C3H7CH=CH2 + CH2=CHOH
    CASE(171)
      j(1016) = seval(szabin,theta,tmp,tmp2,b,c,d) ! n-C6H13CHO -> 2-propylcyclobutanol
    CASE(224)
      j(1017) = seval(szabin,theta,tmp,tmp2,b,c,d) ! nALD -> products
    CASE(141)
      j(1018) = seval(szabin,theta,tmp,tmp2,b,c,d) ! i-C3H7CHO -> i-C3H7 + CHO
    CASE(153)
      j(1019) = seval(szabin,theta,tmp,tmp2,b,c,d) ! i-C4H9CHO -> C4H9 + CHO
    CASE(154)
      j(1020) = seval(szabin,theta,tmp,tmp2,b,c,d) ! i-C4H9CHO -> CH3CH=CH2 + CH2=CHOH
    CASE(155)
      j(1021) = seval(szabin,theta,tmp,tmp2,b,c,d) ! sec-C4H9CHO -> C4H9 + CHO
    CASE(156)
      j(1022) = seval(szabin,theta,tmp,tmp2,b,c,d) ! sec-C4H9CHO  -> CH3CH=CHOH + CH2=CH2
    CASE(157)
      j(1023) = seval(szabin,theta,tmp,tmp2,b,c,d) ! t-C4H9CHO -> C4H9 + CHO
    CASE(185)
      j(1024) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C4H9CH(C2H5)CHO -> C7H15 + CHO
    CASE(186)
      j(1025) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C4H9CH(C2H5)CHO -> C7H16 + CO
    CASE(227)
      j(1026) = seval(szabin,theta,tmp,tmp2,b,c,d) ! MeALD -> products
    CASE(187)
      j(1027) = seval(szabin,theta,tmp,tmp2,b,c,d) ! AlkALD -> products
    CASE(158)
      j(1028) = seval(szabin,theta,tmp,tmp2,b,c,d) ! tALD -> products
    CASE( 52)
      j(1029) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH2=CHCHO -> CH2=CH + CHO
    CASE( 53)
      j(1030) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH2=CHCHO -> CH2=CH2 + CO
    CASE( 55)
      j(1031) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH2=CHCHO -> CH2=CHCO + H
    CASE(179)
      j(1032) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3CH=CHCHO -> CH3CH=CH + CHO
    CASE(180)
      j(1033) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3CH=CHCHO -> CH3CH=CH2 + CO
    CASE(181)
      j(1034) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3CH=CHCHO -> CH3CH=CHCO + H
    CASE(182)
      j(1035) = seval(szabin,theta,tmp,tmp2,b,c,d) ! 2-hexenal -> 1-pentenyl radical + CHO
    CASE(183)
      j(1036) = seval(szabin,theta,tmp,tmp2,b,c,d) ! 2-hexenal -> 1-pentene + CO
    CASE(184)
      j(1037) = seval(szabin,theta,tmp,tmp2,b,c,d) ! 2-hexenal -> C3H7CH=CHCO + H
    CASE(212)
      j(1038) = seval(szabin,theta,tmp,tmp2,b,c,d) ! hexadienal -> 1-pentenyl radical + CHO
    CASE(213)
      j(1039) = seval(szabin,theta,tmp,tmp2,b,c,d) ! hexadienal -> 1,3-pentadiene + CO
    CASE(214)
      j(1040) = seval(szabin,theta,tmp,tmp2,b,c,d) ! hexadienal -> CH3CH=CHCH=CHCO + H
    CASE( 56)
      j(1041) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH2=C(CH3)CHO -> CH2=CCH3 + CHO
    CASE( 57)
      j(1042) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH2=C(CH3)CHO -> CH3CH=CH2 + CO
    CASE( 59)
      j(1043) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH2=C(CH3)CHO -> CH2=C(CH3)CO + H
    CASE(201)
      j(1044) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3C(CH3)=CHCHO -> (CH3)2C=CH + CHO
    CASE(202)
      j(1045) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3C(CH3)=CHCHO -> (CH3)2C=CH2 + CO
    CASE(203)
      j(1046) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3C(CH3)=CHCHO -> (CH3)2C=CHCO + H
    CASE(190)
      j(1047) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3CH=C(CH3)CHO -> CH3CH=CCH3 + CHO
    CASE(191)
      j(1048) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3CH=C(CH3)CHO -> CH3CH=CHCH3 + CO
    CASE(192)
      j(1049) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3CH=C(CH3)CHO -> CH3CH=C(CH3)CO + H
    CASE(230)
      j(1050) = seval(szabin,theta,tmp,tmp2,b,c,d) ! luALD -> NI products
    CASE(231)
      j(1051) = seval(szabin,theta,tmp,tmp2,b,c,d) ! luALD -> alkene + CO
    CASE(232)
      j(1052) = seval(szabin,theta,tmp,tmp2,b,c,d) ! luALD -> acyl + H
    CASE( 64)
      j(1053) = seval(szabin,theta,tmp,tmp2,b,c,d) ! HOCH2CHO -> products
    CASE(178)
      j(1054) = seval(szabin,theta,tmp,tmp2,b,c,d) ! Glycidaldehyde -> products
    CASE( 31)
      j(1055) = seval(szabin,theta,tmp,tmp2,b,c,d) ! ALD3OHqy -> R(OH) + HCO
    CASE(137)
      j(1056) = seval(szabin,theta,tmp,tmp2,b,c,d) ! ALD4OHqy -> NI products
    CASE(138)
      j(1057) = seval(szabin,theta,tmp,tmp2,b,c,d) ! ALD4OHqy -> NII products
    CASE(147)
      j(1058) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C5nALDOHqy -> NI products
    CASE(148)
      j(1059) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C5nALDOHqy -> NII products
    CASE(149)
      j(1060) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C5nALDOHqy -> cycl. product
    CASE(163)
      j(1061) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C6nALDOHqy -> NI products
    CASE(164)
      j(1062) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C6nALDOHqy -> NII products
    CASE(165)
      j(1063) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C6nALDOHqy -> cycl. product
    CASE(172)
      j(1064) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C7nALDOHqy -> NI products
    CASE(173)
      j(1065) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C7nALDOHqy -> NII products
    CASE(174)
      j(1066) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C7nALDOHqy -> cycl. product
    CASE(225)
      j(1067) = seval(szabin,theta,tmp,tmp2,b,c,d) ! nALDOHqy -> products
    CASE(228)
      j(1068) = seval(szabin,theta,tmp,tmp2,b,c,d) ! MeALDOHqy -> products
    CASE(188)
      j(1069) = seval(szabin,theta,tmp,tmp2,b,c,d) ! AlkALDOHqy -> products
    CASE(159)
      j(1070) = seval(szabin,theta,tmp,tmp2,b,c,d) ! tALDOH -> products
    CASE(233)
      j(1071) = seval(szabin,theta,tmp,tmp2,b,c,d) ! luALDOHqy -> NI products
    CASE(234)
      j(1072) = seval(szabin,theta,tmp,tmp2,b,c,d) ! luALDOHqy -> alkene + CO
    CASE(235)
      j(1073) = seval(szabin,theta,tmp,tmp2,b,c,d) ! luALDOHqy -> acyl + H
    CASE(215)
      j(1074) = seval(szabin,theta,tmp,tmp2,b,c,d) ! uuALDOHqy -> NI products
    CASE(216)
      j(1075) = seval(szabin,theta,tmp,tmp2,b,c,d) ! uuALDOHqy -> diene + CO
    CASE(217)
      j(1076) = seval(szabin,theta,tmp,tmp2,b,c,d) ! uuALDOHqy -> acyl + H
    CASE(193)
      j(1077) = seval(szabin,theta,tmp,tmp2,b,c,d) ! aMeC4uALDOHqy -> NI products
    CASE(194)
      j(1078) = seval(szabin,theta,tmp,tmp2,b,c,d) ! aMeC4uALDOHqy -> alkene + CO
    CASE(195)
      j(1079) = seval(szabin,theta,tmp,tmp2,b,c,d) ! aMeC4uALDOHqy -> acyl + H
    CASE(204)
      j(1080) = seval(szabin,theta,tmp,tmp2,b,c,d) ! bMeC4uALDOHqy -> NI products
    CASE(205)
      j(1081) = seval(szabin,theta,tmp,tmp2,b,c,d) ! bMeC4uALDOHqy -> alkene + CO
    CASE(206)
      j(1082) = seval(szabin,theta,tmp,tmp2,b,c,d) ! bMeC4uALDOHqy -> acyl + H
    CASE( 65)
      j(2001) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3COCH3 -> CH3CO + CH3
    CASE( 66)
      j(2002) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3COCH3 -> CO + 2 CH3
    CASE( 67)
      j(2003) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3COCH2CH3 -> CH3CO + CH2CH3
    CASE(243)
      j(2004) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C3H7COCH3 -> CH3CO + C3H7
    CASE(245)
      j(2005) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C3H7COCH3 -> C3H7 + CO + CH3
    CASE(246)
      j(2006) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C3H7COCH3 -> CH3C(OH)=CH2 + CH2=CH2
    CASE(242)
      j(2007) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C2H5COC2H5 -> C2H5CO + C2H5
    CASE(241)
      j(2008) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C2H5COC2H5 -> 2 C2H5 + CO
    CASE(247)
      j(2009) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C4H9COCH3 -> CH3CH=CH2 + CH2=C(OH)CH3
    CASE(248)
      j(2010) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C3H7COC2H5 -> C2H5CO + C3H7
    CASE(249)
      j(2011) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C3H7COC2H5 -> C3H7CO + C2H5
    CASE(250)
      j(2012) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C3H7COC2H5 -> C3H7 + CO + C2H5
    CASE(251)
      j(2013) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C3H7COC2H5 -> C2H5C(OH)=CH2 + CH2=CH2
    CASE(293)
      j(2014) = seval(szabin,theta,tmp,tmp2,b,c,d) ! lcKet -> products
    CASE(252)
      j(2015) = seval(szabin,theta,tmp,tmp2,b,c,d) ! MIPK -> CH3CO + i-C3H7
    CASE(253)
      j(2016) = seval(szabin,theta,tmp,tmp2,b,c,d) ! MIPK -> i-C3H7CO + CH3
    CASE(256)
      j(2017) = seval(szabin,theta,tmp,tmp2,b,c,d) ! MIBK -> products
    CASE(264)
      j(2018) = seval(szabin,theta,tmp,tmp2,b,c,d) ! 5-Me-2-hexanone -> products
    CASE(262)
      j(2019) = seval(szabin,theta,tmp,tmp2,b,c,d) ! 4-Me-2-hexanone -> CH3C(OH)=CH2 + 2-butene
    CASE(263)
      j(2020) = seval(szabin,theta,tmp,tmp2,b,c,d) ! 4-Me-2-hexanone -> CH3C(OH)=CH2 + 1-butene
    CASE(269)
      j(2021) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3CH(CH3)COCH(CH3)2 -> products
    CASE(270)
      j(2022) = seval(szabin,theta,tmp,tmp2,b,c,d) ! c-C3H4O -> C2H4 + CO
    CASE(271)
      j(2023) = seval(szabin,theta,tmp,tmp2,b,c,d) ! c-C3H4O -> further products
    CASE(273)
      j(2024) = seval(szabin,theta,tmp,tmp2,b,c,d) ! c-C4H6O -> C2H4 + CH2=C=O
    CASE(275)
      j(2025) = seval(szabin,theta,tmp,tmp2,b,c,d) ! c-C4H6O -> c-C3H6 + CO
    CASE(277)
      j(2026) = seval(szabin,theta,tmp,tmp2,b,c,d) ! c-C5H8O -> 2 C2H4 + CO
    CASE(278)
      j(2027) = seval(szabin,theta,tmp,tmp2,b,c,d) ! c-C5H8O -> c-C4H8 + CO
    CASE(279)
      j(2028) = seval(szabin,theta,tmp,tmp2,b,c,d) ! c-C5H8O -> CH2=CHCH2CH2CHO
    CASE(281)
      j(2029) = seval(szabin,theta,tmp,tmp2,b,c,d) ! c-C6H10O -> 5-hexenal
    CASE(283)
      j(2030) = seval(szabin,theta,tmp,tmp2,b,c,d) ! c-C6H10O -> 1-pentene + CO
    CASE( 61)
      j(2031) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3COCH=CH2 -> CH3 + C2H3CO
    CASE( 62)
      j(2032) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3COCH=CH2 -> C2H3 + CH3CO
    CASE( 63)
      j(2033) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3COCH=CH2 -> C3H6 + CO
    CASE(285)
      j(2034) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3CH2COCH=CH2 -> C2H5 + C2H3CO
    CASE(286)
      j(2035) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3CH2COCH=CH2 -> C2H3 + C2H5CO
    CASE(287)
      j(2036) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3CH2COCH=CH2 -> 1-C4H8 + CO
    CASE(308)
      j(2037) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH2=C=O -> CO2 + CO + H2
    CASE(309)
      j(2038) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3CH=C=O -> C2H4 + CO
    CASE( 68)
      j(2039) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH2(OH)COCH3 -> CH3CO + CH2(OH)
    CASE( 69)
      j(2040) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH2(OH)COCH3 -> CH2(OH)CO + CH3
    CASE(289)
      j(2041) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3COC2H4OH -> CH3CO + CH2CH2OH
    CASE(290)
      j(2042) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3COCH(OH)CH3 -> CH3CO + CH3CHOH
    CASE(292)
      j(2043) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3COC(CH3)2OH -> CH3CO + (CH3)2COH
    CASE(294)
      j(2044) = seval(szabin,theta,tmp,tmp2,b,c,d) ! lcKet(OH) -> products
    CASE(295)
      j(2045) = seval(szabin,theta,tmp,tmp2,b,c,d) ! cKet(t-OH) -> products
    CASE(297)
      j(2046) = seval(szabin,theta,tmp,tmp2,b,c,d) ! a-br. Ket(OH) -> products
    CASE(298)
      j(2047) = seval(szabin,theta,tmp,tmp2,b,c,d) ! a-br. Ket(t-OH) -> products
    CASE(299)
      j(2048) = seval(szabin,theta,tmp,tmp2,b,c,d) ! b-br. Ket(OH) -> products
    CASE(303)
      j(2049) = seval(szabin,theta,tmp,tmp2,b,c,d) ! uKet(t-OH) -> RdCO + ROH
    CASE(304)
      j(2050) = seval(szabin,theta,tmp,tmp2,b,c,d) ! uKet(t-OH) -> ROHCO + Rd
    CASE(305)
      j(2051) = seval(szabin,theta,tmp,tmp2,b,c,d) ! uKet(t-OH) -> alkene + CO
    CASE(300)
      j(2052) = seval(szabin,theta,tmp,tmp2,b,c,d) ! uKet(OH) -> RdCO + ROH
    CASE(301)
      j(2053) = seval(szabin,theta,tmp,tmp2,b,c,d) ! uKet(OH) -> ROHCO + Rd
    CASE(302)
      j(2054) = seval(szabin,theta,tmp,tmp2,b,c,d) ! uKet(OH) -> alkene + CO
    CASE( 83)
      j(2055) = seval(szabin,theta,tmp,tmp2,b,c,d) ! genPYR(poly)
    CASE( 71)
      j(3001) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CHOCHO -> H2 + 2 CO
    CASE( 72)
      j(3002) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CHOCHO -> CH2O + CO
    CASE( 70)
      j(3003) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CHOCHO -> 2 HO2 + 2 CO
    CASE( 74)
      j(3004) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3COCHO -> CH3CO + HCO
    CASE( 75)
      j(3005) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3COCOCH3 -> Products
    CASE(311)
      j(3006) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CHOCH=CHCHO -> 3H-furan-2-one
    CASE(313)
      j(3007) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3COCH=CHCHO -> 5Me-3H-2-furanone
    CASE(314)
      j(3008) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3COCH=CHCHO -> CH3 + CHOCH=CHCO
    CASE(315)
      j(3009) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3COCH=CHCHO -> CH3COCH=CH2 + CO
    CASE(319)
      j(3010) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3COCH=CHCOCH3 -> CH3CO + CH=CHCOCH3
    CASE(318)
      j(3011) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CHOCH=CHCH=CHCHO -> diformyl cyclobutene
    CASE(142)
      j(3012) = seval(szabin,theta,tmp,tmp2,b,c,d) ! pinonaldehyde -> R + CO + HO2
    CASE(143)
      j(3013) = seval(szabin,theta,tmp,tmp2,b,c,d) ! caronaldehyde -> R + CO + HO2
    CASE( 35)
      j(4001) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3ONO2 -> CH3O + NO2
    CASE( 39)
      j(4002) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C2H5ONO2 -> C2H5O + NO2
    CASE( 40)
      j(4003) = seval(szabin,theta,tmp,tmp2,b,c,d) ! n-C3H7ONO2 -> C3H7O + NO2
    CASE( 41)
      j(4004) = seval(szabin,theta,tmp,tmp2,b,c,d) ! 1-C4H9ONO2 -> 1-C4H9O + NO2
    CASE(320)
      j(4005) = seval(szabin,theta,tmp,tmp2,b,c,d) ! n-C5H11ONO2 -> n-C5H11O + NO2
    CASE( 43)
      j(4006) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3CHONO2CH3 -> CH3CHOCH3 + NO2
    CASE( 42)
      j(4007) = seval(szabin,theta,tmp,tmp2,b,c,d) ! 2-C4H9ONO2 -> 2-C4H9O + NO2
    CASE(321)
      j(4008) = seval(szabin,theta,tmp,tmp2,b,c,d) ! 2-C5H11ONO2 -> 2-C5H11O + NO2
    CASE(322)
      j(4009) = seval(szabin,theta,tmp,tmp2,b,c,d) ! 3-C5H11ONO2 -> 3-C5H11O + NO2
    CASE(323)
      j(4010) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C5H11ONO2 -> C5H11O + NO2
    CASE(325)
      j(4011) = seval(szabin,theta,tmp,tmp2,b,c,d) ! i-C4H9ONO2 -> i-C4H9O + NO2
    CASE( 46)
      j(4012) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C(CH3)3(ONO2) -> C(CH3)3(O.) + NO2
    CASE(326)
      j(4013) = seval(szabin,theta,tmp,tmp2,b,c,d) ! i-C5H11ONO2 -> i-C5H11O + NO2
    CASE(324)
      j(4014) = seval(szabin,theta,tmp,tmp2,b,c,d) ! c-C5H11ONO2 -> c-C5H11O + NO2
    CASE( 44)
      j(4015) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH2(OH)CH2(ONO2) -> CH2(OH)CH2(O.) + NO2
    CASE(328)
      j(4016) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C1(OH)NO3 -> C1(OH)O + NO2
    CASE(329)
      j(4017) = seval(szabin,theta,tmp,tmp2,b,c,d) ! R(OH)NO3 -> R(OH)O + NO2
    CASE(330)
      j(4018) = seval(szabin,theta,tmp,tmp2,b,c,d) ! iR(OH)NO3 -> iR(OH)O + NO2
    CASE(331)
      j(4019) = seval(szabin,theta,tmp,tmp2,b,c,d) ! tR(OH)NO3 -> tR(OH)O + NO2
    CASE(353)
      j(4020) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3NO2 -> CH3 + NO2
    CASE(354)
      j(4021) = seval(szabin,theta,tmp,tmp2,b,c,d) ! RNO2 -> 1-alkene + HONO
    CASE( 38)
      j(4022) = seval(szabin,theta,tmp,tmp2,b,c,d) ! ROONO2 -> products
    CASE( 48)
      j(4023) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3CO(OONO2) -> CH3CO(OO) + NO2
    CASE( 49)
      j(4024) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3CO(OONO2) -> CH3CO(O) + NO3
    CASE( 50)
      j(4025) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3CH2CO(OONO2) -> CH3CH2CO(OO) + NO2
    CASE( 51)
      j(4026) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3CH2CO(OONO2) -> CH3CH2CO(O) + NO3
    CASE(332)
      j(4027) = seval(szabin,theta,tmp,tmp2,b,c,d) ! PAN -> RCO(OO) + NO2
    CASE(333)
      j(4028) = seval(szabin,theta,tmp,tmp2,b,c,d) ! PAN -> RCO(O) + NO3
    CASE(334)
      j(5001) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3CH(NO3)CH2NO3 -> products
    CASE(335)
      j(5002) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C2H5CH(NO3)CH2NO3 -> products
    CASE(336)
      j(5003) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3CH(NO3)CH(NO3)CH3 -> RO. + NO2
    CASE(340)
      j(5004) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C6H9-1-CH3-1,2-NO3 -> products
    CASE(337)
      j(5005) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH2(NO3)CH=CHCH2NO3 -> RO. + NO2
    CASE(338)
      j(5006) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH2=CHCH(NO3)CH2NO3 -> products
    CASE(339)
      j(5007) = seval(szabin,theta,tmp,tmp2,b,c,d) ! uDINIT -> RO. + NO2
    CASE( 33)
      j(6001) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3OOH -> CH3O + OH
    CASE(341)
      j(6002) = seval(szabin,theta,tmp,tmp2,b,c,d) ! (CH3)3COOH -> (CH3)3CO + OH
    CASE( 34)
      j(6003) = seval(szabin,theta,tmp,tmp2,b,c,d) ! HOCH2OOH -> HOCH2O. + OH
    CASE( 78)
      j(6004) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3CO(OOH) -> Products
    CASE(347)
      j(7001) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH2OO -> HCHO + O(3P)
    CASE(348)
      j(7002) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3CHOO -> CH3CHO + O(3P)
    CASE(349)
      j(7003) = seval(szabin,theta,tmp,tmp2,b,c,d) ! synCH3CHOO -> CH3CHO + O(3P)
    CASE(350)
      j(7004) = seval(szabin,theta,tmp,tmp2,b,c,d) ! antiCH3CHOO -> CH3CHO + O(3P)
    CASE(351)
      j(7005) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C2H5CHOO -> C2H5CHO + O(3P)
    CASE(352)
      j(7006) = seval(szabin,theta,tmp,tmp2,b,c,d) ! (CH3)2COO -> CH3COCH3 + O(3P)
    CASE( 45)
      j(8001) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3COCH2(ONO2) -> CH3COCH2(O.) + NO2
    CASE(342)
      j(8002) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C2H5COCH2NO3 -> C2H5COCH2O + NO2
    CASE(343)
      j(8003) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3COCH(NO3)CH3 -> CH3COCH(O.)CH3 + NO2
    CASE(344)
      j(8004) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C6H9-2-=O-1-NO3 -> RO. + NO2
    CASE(345)
      j(8005) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3COCH2CH2CH(OOH)CH3 -> RO. + OH
    CASE(346)
      j(8006) = seval(szabin,theta,tmp,tmp2,b,c,d) ! oxohexyl-hydroperoxide -> RO. + OH
    CASE( 27)
      j(8007) = seval(szabin,theta,tmp,tmp2,b,c,d) ! genCH3CHO(poly)
    CASE( 29)
      j(8008) = seval(szabin,theta,tmp,tmp2,b,c,d) ! genC2H5CHO(poly)
    CASE(239)
      j(8009) = seval(szabin,theta,tmp,tmp2,b,c,d) ! genluALD(poly)
    CASE(210)
      j(8010) = seval(szabin,theta,tmp,tmp2,b,c,d) ! genbMeuAld(poly)
    CASE(199)
      j(8011) = seval(szabin,theta,tmp,tmp2,b,c,d) ! genaMeuAld(poly)
    CASE(221)
      j(8012) = seval(szabin,theta,tmp,tmp2,b,c,d) ! genluuALD(poly)
    CASE( 30)
      j(8013) = seval(szabin,theta,tmp,tmp2,b,c,d) ! genC2H5CHO(OHpoly)
    CASE(240)
      j(8014) = seval(szabin,theta,tmp,tmp2,b,c,d) ! genluALD(OHpoly)
    CASE(211)
      j(8015) = seval(szabin,theta,tmp,tmp2,b,c,d) ! genbMeuAldOH(poly)
    CASE(200)
      j(8016) = seval(szabin,theta,tmp,tmp2,b,c,d) ! genaMeuAldOH(poly)
    CASE(222)
      j(8017) = seval(szabin,theta,tmp,tmp2,b,c,d) ! genluuALD(OHpoly)
    CASE(272)
      j(8018) = seval(szabin,theta,tmp,tmp2,b,c,d) ! genC3cKet
    CASE(276)
      j(8019) = seval(szabin,theta,tmp,tmp2,b,c,d) ! genC4cKet
    CASE(280)
      j(8020) = seval(szabin,theta,tmp,tmp2,b,c,d) ! genC5cKet
    CASE(284)
      j(8021) = seval(szabin,theta,tmp,tmp2,b,c,d) ! genC6cKet
    CASE(288)
      j(8022) = seval(szabin,theta,tmp,tmp2,b,c,d) ! genuKet(poly)
    CASE(310)
      j(8023) = seval(szabin,theta,tmp,tmp2,b,c,d) ! genKete(poly)
    CASE(306)
      j(8024) = seval(szabin,theta,tmp,tmp2,b,c,d) ! genuKet(OHpoly)
    CASE(307)
      j(8025) = seval(szabin,theta,tmp,tmp2,b,c,d) ! genuKet(t-OHpoly)
    CASE( 73)
      j(8027) = seval(szabin,theta,tmp,tmp2,b,c,d) ! genDICAR(Ald)
    CASE( 76)
      j(8028) = seval(szabin,theta,tmp,tmp2,b,c,d) ! genDICAR(Ket)
  END SELECT
