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
      j( 101) = seval(n,theta,tmp,tmp2,b,c,d) ! CH2O -> H + HCO
    CASE( 23)
      j( 102) = seval(n,theta,tmp,tmp2,b,c,d) ! CH2O -> H2 + CO
    CASE( 24)
      j( 103) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3CHO -> CH3 + HCO
    CASE( 26)
      j( 104) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3CHO -> CH3CO + H
    CASE( 27)
      j( 105) = seval(n,theta,tmp,tmp2,b,c,d) ! C2H5CHO -> C2H5 + HCO
    CASE( 28)
      j( 161) = seval(n,theta,tmp,tmp2,b,c,d) ! ALD3OHqy -> R(OH) + HCO
    CASE( 30)
      j(  41) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3OOH -> CH3O + OH
    CASE( 32)
      j(  51) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3ONO2 -> CH3O + NO2
    CASE( 36)
      j(  52) = seval(n,theta,tmp,tmp2,b,c,d) ! C2H5ONO2 -> C2H5O + NO2
    CASE( 37)
      j(  53) = seval(n,theta,tmp,tmp2,b,c,d) ! n-C3H7ONO2 -> C3H7O + NO2
    CASE( 40)
      j(  54) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3CHONO2CH3 -> CH3CHOCH3 + NO2
    CASE( 42)
      j(  56) = seval(n,theta,tmp,tmp2,b,c,d)*0.750 ! CH3COCH2(ONO2) -> CH3COCH2(O.) + NO2
    CASE( 43)
      j(  55) = seval(n,theta,tmp,tmp2,b,c,d) ! C(CH3)3(ONO2) -> C(CH3)3(O.) + NO2
    CASE( 49)
      j( 132) = seval(n,theta,tmp,tmp2,b,c,d) ! CH2=CHCHO -> CH2=CH + CHO
    CASE( 50)
      j( 133) = seval(n,theta,tmp,tmp2,b,c,d) ! CH2=CHCHO -> CH2=CH2 + CO
    CASE( 52)
      j( 134) = seval(n,theta,tmp,tmp2,b,c,d) ! CH2=CHCHO -> CH2=CHCO + H
    CASE( 53)
      j( 144) = seval(n,theta,tmp,tmp2,b,c,d) ! CH2=C(CH3)CHO -> CH2=CCH3 + CHO
    CASE( 54)
      j( 145) = seval(n,theta,tmp,tmp2,b,c,d) ! CH2=C(CH3)CHO -> CH3CH=CH2 + CO
    CASE( 56)
      j( 146) = seval(n,theta,tmp,tmp2,b,c,d) ! CH2=C(CH3)CHO -> CH2=C(CH3)CO + H
    CASE( 59)
      j(  24) = seval(n,theta,tmp,tmp2,b,c,d)*0.500 ! CH3COCH=CH2 -> C2H3 + CH3CO
      j(  23) = seval(n,theta,tmp,tmp2,b,c,d)*0.500 ! CH3COCH=CH2 -> C2H3 + CH3CO
    CASE( 61)
      j( 156) = seval(n,theta,tmp,tmp2,b,c,d) ! HOCH2CHO -> CH2OH + HCO
    CASE( 62)
      j( 157) = seval(n,theta,tmp,tmp2,b,c,d) ! HOCH2CHO -> CH3OH + CO
    CASE( 63)
      j( 158) = seval(n,theta,tmp,tmp2,b,c,d) ! HOCH2CHO -> CH2CHO + OH
    CASE( 64)
      j(  21) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3COCH3 -> CH3CO + CH3
    CASE( 66)
      j(  22) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3COCH2CH3 -> CH3CO + CH2CH3
    CASE( 69)
      j(  33) = seval(n,theta,tmp,tmp2,b,c,d) ! CHOCHO -> 2 HO2 + 2 CO
    CASE( 70)
      j(  31) = seval(n,theta,tmp,tmp2,b,c,d) ! CHOCHO -> H2 + 2 CO
    CASE( 71)
      j(  32) = seval(n,theta,tmp,tmp2,b,c,d) ! CHOCHO -> CH2O + CO
    CASE( 72)
      j(  34) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3COCHO -> CH3CO + HCO
    CASE( 73)
      j(  35) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3COCOCH3 -> Products
    CASE( 82)
      j(1008) = seval(n,theta,tmp,tmp2,b,c,d) ! Cl2 -> Cl + Cl
    CASE( 93)
      j(1006) = seval(n,theta,tmp,tmp2,b,c,d) ! ClONO2 -> Cl + NO3
    CASE( 94)
      j(1007) = seval(n,theta,tmp,tmp2,b,c,d) ! ClONO2 -> ClO + NO2
    CASE(114)
      j(1003) = seval(n,theta,tmp,tmp2,b,c,d) ! Br2 -> Br + Br
    CASE(115)
      j(1002) = seval(n,theta,tmp,tmp2,b,c,d) ! BrO -> Br + O
    CASE(116)
      j(1001) = seval(n,theta,tmp,tmp2,b,c,d) ! HOBr -> OH + Br
    CASE(121)
      j(1005) = seval(n,theta,tmp,tmp2,b,c,d) ! BrONO2 -> BrO + NO2
    CASE(122)
      j(1004) = seval(n,theta,tmp,tmp2,b,c,d) ! BrONO2 -> Br + NO3
    CASE(131)
      j( 106) = seval(n,theta,tmp,tmp2,b,c,d) ! n-C3H7CHO -> n-C3H7 + CHO
    CASE(132)
      j( 107) = seval(n,theta,tmp,tmp2,b,c,d) ! n-C3H7CHO -> C2H4 + CH2CHOH
    CASE(133)
      j( 162) = seval(n,theta,tmp,tmp2,b,c,d) ! ALD4OHqy -> NI products
    CASE(134)
      j( 163) = seval(n,theta,tmp,tmp2,b,c,d) ! ALD4OHqy -> NII products
    CASE(137)
      j( 118) = seval(n,theta,tmp,tmp2,b,c,d) ! i-C3H7CHO -> i-C3H7 + CHO
    CASE(140)
      j( 108) = seval(n,theta,tmp,tmp2,b,c,d) ! n-C4H9CHO -> C4H9 +  CHO
    CASE(141)
      j( 109) = seval(n,theta,tmp,tmp2,b,c,d) ! n-C4H9CHO -> CH3CH=CH2 +  CH2=CHOH
    CASE(142)
      j( 110) = seval(n,theta,tmp,tmp2,b,c,d) ! n-C4H9CHO -> 2-methylcyclobutanol
    CASE(143)
      j( 164) = seval(n,theta,tmp,tmp2,b,c,d) ! C5nALDOHqy -> NI products
    CASE(144)
      j( 165) = seval(n,theta,tmp,tmp2,b,c,d) ! C5nALDOHqy -> NII products
    CASE(148)
      j( 119) = seval(n,theta,tmp,tmp2,b,c,d) ! i-C4H9CHO -> C4H9 + CHO
    CASE(149)
      j( 120) = seval(n,theta,tmp,tmp2,b,c,d) ! i-C4H9CHO -> CH3CH=CH2 + CH2=CHOH
    CASE(150)
      j( 121) = seval(n,theta,tmp,tmp2,b,c,d) ! sec-C4H9CHO -> C4H9 + CHO
    CASE(151)
      j( 122) = seval(n,theta,tmp,tmp2,b,c,d) ! sec-C4H9CHO  -> CH3CH=CHOH + CH2=CH2
    CASE(152)
      j( 123) = seval(n,theta,tmp,tmp2,b,c,d) ! t-C4H9CHO -> C4H9 + CHO
    CASE(153)
      j( 130) = seval(n,theta,tmp,tmp2,b,c,d) ! tALD -> NI products
    CASE(154)
      j( 131) = seval(n,theta,tmp,tmp2,b,c,d) ! tALD -> NII products
    CASE(155)
      j( 175) = seval(n,theta,tmp,tmp2,b,c,d) ! tALDOHqy -> NI products
    CASE(156)
      j( 176) = seval(n,theta,tmp,tmp2,b,c,d) ! tALDOHqy -> NII products
    CASE(159)
      j( 111) = seval(n,theta,tmp,tmp2,b,c,d) ! n-C5H11CHO -> C5H11 + CHO
    CASE(160)
      j( 112) = seval(n,theta,tmp,tmp2,b,c,d) ! n-C5H11CHO -> C2H5CH=CH2 + CH2=CHOH
    CASE(161)
      j( 113) = seval(n,theta,tmp,tmp2,b,c,d) ! n-C5H11CHO -> 2-ethylcyclobutanol
    CASE(162)
      j( 166) = seval(n,theta,tmp,tmp2,b,c,d) ! C6nALDOHqy -> NI products
    CASE(163)
      j( 167) = seval(n,theta,tmp,tmp2,b,c,d) ! C6nALDOHqy -> NII products
    CASE(167)
      j( 114) = seval(n,theta,tmp,tmp2,b,c,d) ! n-C6H13CHO -> C6H13 + CHO
    CASE(168)
      j( 115) = seval(n,theta,tmp,tmp2,b,c,d) ! n-C6H13CHO -> C3H7CH=CH2 + CH2=CHOH
    CASE(169)
      j( 116) = seval(n,theta,tmp,tmp2,b,c,d) ! n-C6H13CHO -> 2-propylcyclobutanol
    CASE(170)
      j( 168) = seval(n,theta,tmp,tmp2,b,c,d) ! C7nALDOHqy -> NI products
    CASE(171)
      j( 169) = seval(n,theta,tmp,tmp2,b,c,d) ! C7nALDOHqy -> NII products
    CASE(175)
      j( 159) = seval(n,theta,tmp,tmp2,b,c,d) ! Glycidaldehyde -> oxyranyl radical + CHO
    CASE(176)
      j( 160) = seval(n,theta,tmp,tmp2,b,c,d) ! Glycidaldehyde -> oxyrane + CO
    CASE(177)
      j( 135) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3CH=CHCHO -> CH3CH=CH + CHO
    CASE(178)
      j( 136) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3CH=CHCHO -> CH3CH=CH2 + CO
    CASE(179)
      j( 137) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3CH=CHCHO -> CH3CH=CHCO + H
    CASE(180)
      j( 138) = seval(n,theta,tmp,tmp2,b,c,d) ! 2-hexenal -> 1-pentenyl radical + CHO
    CASE(181)
      j( 139) = seval(n,theta,tmp,tmp2,b,c,d) ! 2-hexenal -> 1-pentene + CO
    CASE(182)
      j( 140) = seval(n,theta,tmp,tmp2,b,c,d) ! 2-hexenal -> C3H7CH=CHCO + H
    CASE(183)
      j( 124) = seval(n,theta,tmp,tmp2,b,c,d) ! C4H9CH(C2H5)CHO -> C7H15 + CHO
    CASE(184)
      j( 125) = seval(n,theta,tmp,tmp2,b,c,d) ! C4H9CH(C2H5)CHO -> C7H16 + CO
    CASE(185)
      j( 128) = seval(n,theta,tmp,tmp2,b,c,d) ! AlkALD -> NI products
    CASE(186)
      j( 129) = seval(n,theta,tmp,tmp2,b,c,d) ! AlkALD -> NII products
    CASE(187)
      j( 173) = seval(n,theta,tmp,tmp2,b,c,d) ! AlkALDOHqy -> NI products
    CASE(188)
      j( 174) = seval(n,theta,tmp,tmp2,b,c,d) ! AlkALDOHqy -> NII products
    CASE(191)
      j( 150) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3CH=C(CH3)CHO -> CH3CH=CCH3 + CHO
    CASE(192)
      j( 151) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3CH=C(CH3)CHO -> CH3CH=CHCH3 + CO
    CASE(193)
      j( 152) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3CH=C(CH3)CHO -> CH3CH=C(CH3)CO + H
    CASE(194)
      j( 183) = seval(n,theta,tmp,tmp2,b,c,d) ! aMeC4uALDOHqy -> NI products
    CASE(195)
      j( 184) = seval(n,theta,tmp,tmp2,b,c,d) ! aMeC4uALDOHqy -> alkene + CO
    CASE(196)
      j( 185) = seval(n,theta,tmp,tmp2,b,c,d) ! aMeC4uALDOHqy -> acyl + H
    CASE(200)
      j( 147) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3C(CH3)=CHCHO -> (CH3)2C=CH + CHO
    CASE(201)
      j( 148) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3C(CH3)=CHCHO -> (CH3)2C=CH2 + CO
    CASE(202)
      j( 149) = seval(n,theta,tmp,tmp2,b,c,d) ! CH3C(CH3)=CHCHO -> (CH3)2C=CHCO + H
    CASE(203)
      j( 186) = seval(n,theta,tmp,tmp2,b,c,d) ! bMeC4uALDOHqy -> NI products
    CASE(204)
      j( 187) = seval(n,theta,tmp,tmp2,b,c,d) ! bMeC4uALDOHqy -> alkene + CO
    CASE(205)
      j( 188) = seval(n,theta,tmp,tmp2,b,c,d) ! bMeC4uALDOHqy -> acyl + H
    CASE(209)
      j( 141) = seval(n,theta,tmp,tmp2,b,c,d) ! hexadienal -> 1-pentenyl radical + CHO
    CASE(210)
      j( 142) = seval(n,theta,tmp,tmp2,b,c,d) ! hexadienal -> 1,3-pentadiene + CO
    CASE(211)
      j( 143) = seval(n,theta,tmp,tmp2,b,c,d) ! hexadienal -> CH3CH=CHCH=CHCO + H
    CASE(212)
      j( 177) = seval(n,theta,tmp,tmp2,b,c,d) ! uuALDOHqy -> NI products
    CASE(213)
      j( 178) = seval(n,theta,tmp,tmp2,b,c,d) ! uuALDOHqy -> diene + CO
    CASE(214)
      j( 179) = seval(n,theta,tmp,tmp2,b,c,d) ! uuALDOHqy -> acyl + H
    CASE(219)
      j( 117) = seval(n,theta,tmp,tmp2,b,c,d) ! nALD -> products
    CASE(220)
      j( 170) = seval(n,theta,tmp,tmp2,b,c,d) ! nALDOHqy -> products
    CASE(222)
      j( 126) = seval(n,theta,tmp,tmp2,b,c,d) ! MeALD -> NI products
    CASE(223)
      j( 127) = seval(n,theta,tmp,tmp2,b,c,d) ! MeALD -> NII products
    CASE(224)
      j( 171) = seval(n,theta,tmp,tmp2,b,c,d) ! MeALDOHqy -> NI products
    CASE(225)
      j( 172) = seval(n,theta,tmp,tmp2,b,c,d) ! MeALDOHqy -> NII products
    CASE(228)
      j( 189) = seval(n,theta,tmp,tmp2,b,c,d) ! cALD -> NI products
    CASE(229)
      j( 190) = seval(n,theta,tmp,tmp2,b,c,d) ! cALDOHqy -> NI products
    CASE(231)
      j( 153) = seval(n,theta,tmp,tmp2,b,c,d) ! luALD -> NI products
    CASE(232)
      j( 154) = seval(n,theta,tmp,tmp2,b,c,d) ! luALD -> alkene + CO
    CASE(233)
      j( 155) = seval(n,theta,tmp,tmp2,b,c,d) ! luALD -> acyl + H
    CASE(212)
      j( 177) = seval(n,theta,tmp,tmp2,b,c,d) ! uuALDOHqy -> NI products
    CASE(213)
      j( 178) = seval(n,theta,tmp,tmp2,b,c,d) ! uuALDOHqy -> diene + CO
    CASE(214)
      j( 179) = seval(n,theta,tmp,tmp2,b,c,d) ! uuALDOHqy -> acyl + H
  END SELECT
