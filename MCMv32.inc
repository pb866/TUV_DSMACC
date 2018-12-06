  SELECT CASE (jl)
    CASE(2) ! O3 -> O2 + O(1D)
      j(1) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(3) ! O3 -> O2 + O(3P)
      j(2) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(5) ! H2O2 -> 2 OH
      j(3) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(6) ! NO2 -> NO + O(3P)
      j(4) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(7) ! NO3 -> NO + O2
      j(5) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(8) ! NO3 -> NO2 + O(3P)
      j(6) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(12) ! HNO2 -> OH + NO
      j(7) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(13) ! HNO3 -> OH + NO2
      j(8) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(22) ! CH2O -> H + HCO
      j(11) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(23) ! CH2O -> H2 + CO
      j(12) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(24) ! CH3CHO -> CH3 + HCO
      j(13) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(28) ! C2H5CHO -> C2H5 + HCO
      j(14) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(134) ! n-C3H7CHO -> n-C3H7 + CHO
      j(15) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(135) ! n-C3H7CHO -> C2H4 + CH2CHOH
      j(16) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(140) ! i-C3H7CHO -> i-C3H7 + CHO
      j(17) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(54) ! CH2=C(CH3)CHO -> CH2=CCH3 + CHO
      j(18) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(57) ! CH2=C(CH3)CHO -> CH2=C(CH3)CO + H
      j(19) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(58) ! HPALD -> Products
      j(20) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(65) ! CH3COCH3 -> CH3CO + CH3
      j(21) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(67) ! CH3COCH2CH3 -> CH3CO + CH2CH3
      j(22) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(61) ! CH3COCH=CH2 -> C3H6 + CO
      j(23) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(60) ! CH3COCH=CH2 -> C2H3 + CH3CO
      j(24) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(71) ! CHOCHO -> H2 + 2 CO
      j(31) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(72) ! CHOCHO -> CH2O + CO
      j(32) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(70) ! CHOCHO -> 2 HO2 + 2 CO
      j(33) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(75) ! CH3COCHO -> CH3CO + HCO
      j(34) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(76) ! CH3COCOCH3 -> Products
      j(35) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(32) ! CH3OOH -> CH3O + OH
      j(41) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(34) ! CH3ONO2 -> CH3O + NO2
      j(51) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(37) ! C2H5ONO2 -> C2H5O + NO2
      j(52) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(38) ! n-C3H7ONO2 -> C3H7O + NO2
      j(53) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(41) ! CH3CHONO2CH3 -> CH3CHOCH3 + NO2
      j(54) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(44) ! C(CH3)3(ONO2) -> C(CH3)3(O.) + NO2
      j(55) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(43) ! CH3COCH2(ONO2) -> CH3COCH2(O.) + NO2
      j(56) = seval(szabin,theta,tmp,tmp2,b,c,d)*0.750
      j(57) = seval(szabin,theta,tmp,tmp2,b,c,d)*0.250
    CASE(18) ! HNO4 -> HO2 + NO2
      j(1300) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(119) ! HOBr -> OH + Br
      j(1001) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(118) ! BrO -> Br + O
      j(1002) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(117) ! Br2 -> Br + Br
      j(1003) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(125) ! BrONO2 -> Br + NO3
      j(1004) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(124) ! BrONO2 -> BrO + NO2
      j(1005) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(96) ! ClONO2 -> Cl + NO3
      j(1006) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(97) ! ClONO2 -> ClO + NO2
      j(1007) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(85) ! Cl2 -> Cl + Cl
      j(1008) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(94) ! ClNO2 -> Cl + NO2
      j(1009) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(86) ! ClO -> Cl + O(1D)
      j(1010) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(87) ! ClO -> Cl + O(3P)
      j(1011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(91) ! HCl -> H + Cl
      j(1012) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(123) ! BrNO2 -> Br + NO2
      j(1013) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(365) ! IO -> I + O(3P)
      j(1100) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(366) ! HOI -> I + OH
      j(1101) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(367) ! OIO -> I + O2
      j(1102) = seval(szabin,theta,tmp,tmp2,b,c,d)
  END SELECT
