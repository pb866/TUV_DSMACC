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
    CASE( 22)
      j(  11) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH2O -> H + HCO
    CASE( 23)
      j(  12) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH2O -> H2 + CO
    CASE( 24)
      j(  13) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3CHO -> CH3 + HCO
    CASE( 28)
      j(  14) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C2H5CHO -> C2H5 + HCO
    CASE(134)
      j(  15) = seval(szabin,theta,tmp,tmp2,b,c,d) ! n-C3H7CHO -> n-C3H7 + CHO
    CASE(135)
      j(  16) = seval(szabin,theta,tmp,tmp2,b,c,d) ! n-C3H7CHO -> C2H4 + CH2CHOH
    CASE(140)
      j(  17) = seval(szabin,theta,tmp,tmp2,b,c,d) ! i-C3H7CHO -> i-C3H7 + CHO
    CASE( 54)
      j(  18) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH2=C(CH3)CHO -> CH2=CCH3 + CHO
    CASE( 57)
      j(  19) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH2=C(CH3)CHO -> CH2=C(CH3)CO + H
    CASE( 58)
      j(  20) = seval(szabin,theta,tmp,tmp2,b,c,d) ! HPALD -> Products
    CASE( 65)
      j(  21) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3COCH3 -> CH3CO + CH3
    CASE( 67)
      j(  22) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3COCH2CH3 -> CH3CO + CH2CH3
    CASE( 61)
      j(  23) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3COCH=CH2 -> C3H6 + CO
    CASE( 60)
      j(  24) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3COCH=CH2 -> C2H3 + CH3CO
    CASE( 71)
      j(  31) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CHOCHO -> H2 + 2 CO
    CASE( 72)
      j(  32) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CHOCHO -> CH2O + CO
    CASE( 70)
      j(  33) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CHOCHO -> 2 HO2 + 2 CO
    CASE( 75)
      j(  34) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3COCHO -> CH3CO + HCO
    CASE( 76)
      j(  35) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3COCOCH3 -> Products
    CASE( 32)
      j(  41) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3OOH -> CH3O + OH
    CASE( 34)
      j(  51) = seval(szabin,theta,tmp,tmp2,b,c,d) ! CH3ONO2 -> CH3O + NO2
    CASE( 37)
      j(  52) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C2H5ONO2 -> C2H5O + NO2
      j(  54) = seval(szabin,theta,tmp,tmp2,b,c,d) ! C2H5ONO2 -> C2H5O + NO2
    CASE( 38)
      j(  53) = seval(szabin,theta,tmp,tmp2,b,c,d) ! n-C3H7ONO2 -> C3H7O + NO2
    CASE( 40)
      j(  55) = seval(szabin,theta,tmp,tmp2,b,c,d) ! 2-C4H9ONO2 -> 2-C4H9O + NO2
    CASE( 43)
      j(  56) = seval(szabin,theta,tmp,tmp2,b,c,d)*0.750 ! CH3COCH2(ONO2) -> CH3COCH2(O.) + NO2
      j(  57) = seval(szabin,theta,tmp,tmp2,b,c,d)*0.250 ! CH3COCH2(ONO2) -> CH3COCH2(O.) + NO2
! Amendment for Edwards et al. (2014) mechanism
    CASE( 18)
      j(1300)=seval(szabin,theta,tmp, tmp2, b,c,d)    ! HNO4 -> HO2 + NO2
    CASE(119)
      j(1001)=0.!seval(szabin,theta,tmp, tmp2, b,c,d) ! HOBr
    CASE(118)
      j(1002)=0.!seval(szabin,theta,tmp, tmp2, b,c,d) ! BrO
    CASE(117)
      j(1003)=0.!seval(szabin,theta,tmp, tmp2, b,c,d) ! Br2
    CASE(125)
      j(1004)=0.!seval(szabin,theta,tmp, tmp2, b,c,d) ! BrNO3->Br+NO3
    CASE(124)
      j(1005)=0.!seval(szabin,theta,tmp, tmp2, b,c,d) ! BrNO3->BrO+NO2
    CASE( 96)
      j(1006)=0.!seval(szabin,theta,tmp, tmp2, b,c,d) ! ClNO3->Cl+NO3
    CASE( 97)
      j(1007)=0.!seval(szabin,theta,tmp, tmp2, b,c,d) ! ClNO3->ClO+NO2
    CASE( 88)
      j(1008)=seval(szabin,theta,tmp, tmp2, b,c,d) ! Cl2->2Cl
  END SELECT
