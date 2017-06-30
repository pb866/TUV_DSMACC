
SELECT CASE(jl)

    CASE(  1) !O2 -> O + O
        j(  1) =  seval(szabin,theta,tmp,tmp2,b,c,d)!O3 -> O2 + O(1D)
        
    CASE(  2) !O3 -> O + O2
        j(  3) =  seval(szabin,theta,tmp,tmp2,b,c,d)!HO2 -> OH + O
        
    CASE(  3) !O3 -> O(1D) + O2
        j(  2) =  seval(szabin,theta,tmp,tmp2,b,c,d)!O3 -> O2 + O(3P)
        
    CASE(  5) !HO2 -> O + OH
        j(  4) =  seval(szabin,theta,tmp,tmp2,b,c,d)!H2O2 -> 2 OH
        
    CASE(  7) !CH2O -> CO + H + HO2
        j( 22) =  seval(szabin,theta,tmp,tmp2,b,c,d)!CH2O -> H2 + CO
        
    CASE(  8) !CH2O -> CO + H2
        j( 23) =  seval(szabin,theta,tmp,tmp2,b,c,d)!CH3CHO -> CH3 + HCO
        
    CASE(  9) !H2O2 -> OH + OH
        j(  5) =  seval(szabin,theta,tmp,tmp2,b,c,d)!NO2 -> NO + O(3P)
        
    CASE( 10) !MP -> CH2O + HO2 + OH
        j( 33) =  seval(szabin,theta,tmp,tmp2,b,c,d)!HOCH2OOH -> HOCH2O. + OH
        
    CASE( 11) !NO2 -> NO + O
        j(  6) =  seval(szabin,theta,tmp,tmp2,b,c,d)!NO3 -> NO + O2
        
    CASE( 12) !NO3 -> NO2 + O
        j(  8) =  seval(szabin,theta,tmp,tmp2,b,c,d)!N2O -> N2 + O(1D)
        
    CASE( 13) !NO3 -> NO + O2
        j(  7) =  seval(szabin,theta,tmp,tmp2,b,c,d)!NO3 -> NO2 + O(3P)
        
    CASE( 14) !N2O5 -> NO2 + NO3
        j( 11) =  seval(szabin,theta,tmp,tmp2,b,c,d)!HNO2 -> OH + NO
        
    CASE( 15) !HNO2 -> NO + OH
        j( 12) =  seval(szabin,theta,tmp,tmp2,b,c,d)!HNO3 -> OH + NO2
        
    CASE( 16) !HNO3 -> NO2 + OH
        j( 13) =  seval(szabin,theta,tmp,tmp2,b,c,d)!HNO3 -> OH + NO + O(3P)
        
    CASE( 17) !HNO4 -> NO3 + OH
        j( 19) =  seval(szabin,theta,tmp,tmp2,b,c,d)!NO3-(aq) -> NO2(aq) + O-
        
    CASE( 18) !HNO4 -> HO2 + NO2
        j( 18) =  seval(szabin,theta,tmp,tmp2,b,c,d)!HNO4 -> OH + NO3
        
    CASE( 19) !ClNO3 -> Cl + NO3
        j(100) =  seval(szabin,theta,tmp,tmp2,b,c,d)!ClONO2 -> ClO + NO2
        
    CASE( 20) !ClNO3 -> ClO + NO2
        j(101) =  seval(szabin,theta,tmp,tmp2,b,c,d)!CCl4 -> Products
        
    CASE( 21) !ClNO2 -> Cl + NO2
        j( 98) =  seval(szabin,theta,tmp,tmp2,b,c,d)!ClONO -> Cl + NO2
        
    CASE( 22) !Cl2 -> Cl + Cl
        j( 89) =  seval(szabin,theta,tmp,tmp2,b,c,d)!ClO -> Cl + O(1D)
        
    CASE( 23) !Br2 -> Br + Br
        j(121) =  seval(szabin,theta,tmp,tmp2,b,c,d)!BrO -> Br + O
        
    CASE( 26) !Cl2O2 -> Cl + Cl + O2
        j( 94) =  seval(szabin,theta,tmp,tmp2,b,c,d)!HCl -> H + Cl
        
    CASE( 27) !ClO -> Cl + O
        j( 91) =  seval(szabin,theta,tmp,tmp2,b,c,d)!ClOO -> Products
        
    CASE( 28) !BrO -> Br + O
        j(122) =  seval(szabin,theta,tmp,tmp2,b,c,d)!HOBr -> OH + Br
        
    CASE( 29) !BrNO3 -> Br + NO3
        j(129) =  seval(szabin,theta,tmp,tmp2,b,c,d)!BrCl -> Br + Cl
        
    CASE( 30) !BrNO3 -> BrO + NO2
        j(128) =  seval(szabin,theta,tmp,tmp2,b,c,d)!BrONO2 -> Br + NO3
        
    CASE( 31) !BrNO2 -> Br + NO2
        j(127) =  seval(szabin,theta,tmp,tmp2,b,c,d)!BrONO2 -> BrO + NO2
        
    CASE( 32) !HOBr -> Br + OH
        j(123) =  seval(szabin,theta,tmp,tmp2,b,c,d)!BrNO -> Br + NO
        
    CASE( 33) !BrCl -> Br + Cl
        j(130) =  seval(szabin,theta,tmp,tmp2,b,c,d)!CH3Br -> Products
        
    CASE( 36) !N2O -> N2 + O
        j(  9) =  seval(szabin,theta,tmp,tmp2,b,c,d)!N2O5 -> NO3 + NO + O(3P)
        
    CASE( 37) !CFC11 -> Cl + PRODUCTS
        j(119) =  seval(szabin,theta,tmp,tmp2,b,c,d)!CF2ClCF2CHFCl (HCFC-225cb) -> Products
        
    CASE( 38) !CFC12 -> Cl + PRODUCTS
        j(110) =  seval(szabin,theta,tmp,tmp2,b,c,d)!CF2ClCFCl2 (CFC-113) -> Products
        
    CASE( 39) !CFC113 -> Cl + PRODUCTS
        j(111) =  seval(szabin,theta,tmp,tmp2,b,c,d)!CF2ClCF2Cl (CFC-114) -> Products
        
    CASE( 40) !CFC114 -> Cl + PRODUCTS
        j(112) =  seval(szabin,theta,tmp,tmp2,b,c,d)!CF3CF2Cl (CFC-115) -> Products
        
    CASE( 41) !CFC115 -> Cl + PRODUCTS
        j(113) =  seval(szabin,theta,tmp,tmp2,b,c,d)!CHClF2 (HCFC-22) -> Products
        
    CASE( 42) !CCl4 -> Cl + PRODUCTS
        j(102) =  seval(szabin,theta,tmp,tmp2,b,c,d)!CH3OCl -> CH3O + Cl
        
    CASE( 43) !CH3Cl -> Cl + PRODUCTS
        j(105) =  seval(szabin,theta,tmp,tmp2,b,c,d)!CH3CCl3 -> Products
        
    CASE( 44) !CH3CCl3 -> Cl + PRODUCTS
        j(106) =  seval(szabin,theta,tmp,tmp2,b,c,d)!CCl2O -> Products
        
    CASE( 50) !CH3Br -> Br + PRODUCTS
        j(131) =  seval(szabin,theta,tmp,tmp2,b,c,d)!CHBr3 -> Products
        
    CASE( 51) !H1211 -> Br + PRODUCTS
        j(132) =  seval(szabin,theta,tmp,tmp2,b,c,d)!CF2Br2 (Halon-1202) -> Products
        
    CASE( 56) !CHBr3 -> Br + PRODUCTS
        j(132) =  seval(szabin,theta,tmp,tmp2,b,c,d)!CF2Br2 (Halon-1202) -> Products
        
    CASE( 62) !ALD2 -> CH4 + CO
        j( 24) = 0.01 *  seval(szabin,theta,tmp,tmp2,b,c,d)!CH3CHO -> CH4 + CO
        
    CASE( 63) !MVK -> CO + PRPE
        j( 63) =  seval(szabin,theta,tmp,tmp2,b,c,d)!HOCH2CHO -> CH2OH + HCO
        
    CASE( 64) !MVK -> CH2O + CO + HO2 + MCO3
        j( 62) = 0.5 *  seval(szabin,theta,tmp,tmp2,b,c,d)!CH3COCH=CH2 -> C3H6 + CO
        
    CASE( 65) !MVK -> MO2 + RCO3
        j( 62) = 0.5 *  seval(szabin,theta,tmp,tmp2,b,c,d)!CH3COCH=CH2 -> C3H6 + CO
        
    CASE( 71) !MGLY -> CO + HO2 + MCO3
        j( 76) =  seval(szabin,theta,tmp,tmp2,b,c,d)!aDICARak(poly)
        
    CASE( 72) !GLYX -> CO + CO + HO2 + HO2
        j( 73) =  seval(szabin,theta,tmp,tmp2,b,c,d)!CHOCHO -> H2 + 2 CO
        
    CASE( 73) !GLYX -> CO + CO + H2
        j( 74) =  seval(szabin,theta,tmp,tmp2,b,c,d)!CHOCHO -> CH2O + CO
        
    CASE( 74) !GLYX -> CH2O + CO
        j( 75) =  seval(szabin,theta,tmp,tmp2,b,c,d)!CH3COCHO -> CH3CO + HCO
        
    CASE( 75) !HAC -> HCHO + HO2 + MeCO3
        j( 70) =  seval(szabin,theta,tmp,tmp2,b,c,d)!CH2(OH)COCH3 -> CH3CO + CH2(OH)
        
    CASE( 76) !ACET -> MCO3 + MO2
        j( 68) =  seval(szabin,theta,tmp,tmp2,b,c,d)!CH3COCH3 -> CO + 2 CH3
        
    CASE( 77) !ACET -> CO + MO2 + MO2
        j( 69) =  seval(szabin,theta,tmp,tmp2,b,c,d)!CH3COCH2CH3 -> CH3CO + CH2CH3
        
    CASE( 78) !INPN -> HO2 + NO2 + OH + RCHO
        j( 33) =  seval(szabin,theta,tmp,tmp2,b,c,d)!HOCH2OOH -> HOCH2O. + OH
        
    CASE( 79) !PRPN -> HO2 + NO2 + OH + RCHO
        j( 33) =  seval(szabin,theta,tmp,tmp2,b,c,d)!HOCH2OOH -> HOCH2O. + OH
        
    CASE( 80) !ETP -> ALD2 + HO2 + OH
        j( 33) =  seval(szabin,theta,tmp,tmp2,b,c,d)!HOCH2OOH -> HOCH2O. + OH
        
    CASE( 81) !RA3P -> HO2 + OH + RCHO
        j( 33) =  seval(szabin,theta,tmp,tmp2,b,c,d)!HOCH2OOH -> HOCH2O. + OH
        
    CASE( 82) !RB3P -> ACET + HO2 + OH
        j( 33) =  seval(szabin,theta,tmp,tmp2,b,c,d)!HOCH2OOH -> HOCH2O. + OH
        
    CASE( 83) !R4P -> HO2 + OH + RCHO
        j( 33) =  seval(szabin,theta,tmp,tmp2,b,c,d)!HOCH2OOH -> HOCH2O. + OH
        
    CASE( 84) !PP -> ALD2 + CH2O + HO2 + OH
        j( 41) =  seval(szabin,theta,tmp,tmp2,b,c,d)!2-C4H9ONO2 -> 2-C4H9O + NO2
        
    CASE( 85) !RP -> ALD2 + HO2 + OH
        j( 33) =  seval(szabin,theta,tmp,tmp2,b,c,d)!HOCH2OOH -> HOCH2O. + OH
        
    CASE( 88) !ISNP -> HO2 + NO2 + OH + RCHO
        j( 33) =  seval(szabin,theta,tmp,tmp2,b,c,d)!HOCH2OOH -> HOCH2O. + OH
        
    CASE( 91) !MAOP -> CH2O + MCO3 + OH
        j( 33) =  seval(szabin,theta,tmp,tmp2,b,c,d)!HOCH2OOH -> HOCH2O. + OH
        
    CASE( 99) !MAP -> MO2 + OH
        j( 33) =  seval(szabin,theta,tmp,tmp2,b,c,d)!HOCH2OOH -> HOCH2O. + OH
        
    CASE(101) !ClNO2 -> Cl + NO2
        j( 98) =  seval(szabin,theta,tmp,tmp2,b,c,d)!ClONO -> Cl + NO2
        
    CASE(103) !O3 -> O(1D) + O2
        j(  2) =  seval(szabin,theta,tmp,tmp2,b,c,d)!O3 -> O2 + O(3P)
        
END SELECT