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
    CASE(18) ! HNO4 -> HO2 + NO2
      j(9) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(19) ! HNO4 -> OH + NO3
      j(10) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(22) ! CH2O -> H + HCO
      j(11011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(23) ! CH2O -> H2 + CO
      j(11012) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(24) ! CH3CHO -> CH3 + HCO
      j(11021) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(26) ! CH3CHO -> CH3CO + H
      j(11022) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(28) ! C2H5CHO -> C2H5 + HCO
      j(11031) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(134) ! n-C3H7CHO -> n-C3H7 + CHO
      j(11041) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(135) ! n-C3H7CHO -> C2H4 + CH2CHOH
      j(11042) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(141) ! n-C4H9CHO -> C4H9 + CHO
      j(11051) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(159) ! n-C5H11CHO -> C5H11 + CHO
      j(11061) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(167) ! n-C6H13CHO -> C6H13 + CHO
      j(11071) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(212) ! n-C7H15CHO -> C7H15 + CHO
      j(11081) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(140) ! i-C3H7CHO -> i-C3H7 + CHO
      j(12011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(149) ! i-C4H9CHO -> C4H9 + CHO
      j(12021) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(150) ! i-C4H9CHO -> CH3CH=CH2 + CH2=CHOH
      j(12022) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(152) ! sec-C4H9CHO -> C4H9 + CHO
      j(12031) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(153) ! sec-C4H9CHO -> CH3CH=CHOH + CH2=CH2
      j(12032) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(154) ! t-C4H9CHO -> C4H9 + CHO
      j(12041) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(183) ! C4H9CH(C2H5)CHO -> C7H15 + CHO
      j(12051) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(155) ! tALD -> products
      j(12061) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(156) ! neoC5H11CHO -> neoC5H11 + CHO
      j(12071) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(50) ! CH2=CHCHO -> CH2=CH + CHO
      j(13011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(51) ! CH2=CHCHO -> CH2=CH2 + CO
      j(13012) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(53) ! CH2=CHCHO -> CH2=CHCO + H
      j(13013) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(177) ! CH3CH=CHCHO -> CH3CH=CH + CHO
      j(13021) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(178) ! CH3CH=CHCHO -> CH3CH=CH2 + CO
      j(13022) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(179) ! CH3CH=CHCHO -> CH3CH=CHCO + H
      j(13023) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(180) ! 2-hexenal -> 1-pentenyl radical + CHO
      j(13031) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(181) ! 2-hexenal -> 1-pentene + CO
      j(13032) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(182) ! 2-hexenal -> C3H7CH=CHCO + H
      j(13033) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(203) ! hexadienal -> 1-pentenyl radical + CHO
      j(13041) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(204) ! hexadienal -> 1,3-pentadiene + CO
      j(13042) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(205) ! hexadienal -> CH3CH=CHCH=CHCO + H
      j(13043) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(54) ! CH2=C(CH3)CHO -> CH2=CCH3 + CHO
      j(13051) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(55) ! CH2=C(CH3)CHO -> CH3CH=CH2 + CO
      j(13052) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(57) ! CH2=C(CH3)CHO -> CH2=C(CH3)CO + H
      j(13053) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(195) ! CH3C(CH3)=CHCHO -> (CH3)2C=CH + CHO
      j(13061) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(196) ! CH3C(CH3)=CHCHO -> (CH3)2C=CH2 + CO
      j(13062) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(197) ! CH3C(CH3)=CHCHO -> (CH3)2C=CHCO + H
      j(13063) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(187) ! CH3CH=C(CH3)CHO -> CH3CH=CCH3 + CHO
      j(13071) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(188) ! CH3CH=C(CH3)CHO -> CH3CH=CHCH3 + CO
      j(13072) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(189) ! CH3CH=C(CH3)CHO -> CH3CH=C(CH3)CO + H
      j(13073) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(220) ! luALD -> NI products
      j(13081) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(221) ! luALD -> alkene + CO
      j(13082) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(222) ! luALD -> acyl + H
      j(13083) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(62) ! HOCH2CHO -> CH2OH + HCO
      j(15011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(29) ! ALD3OH -> R(OH) + HCO
      j(15021) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(136) ! ALD4OH -> NI products
      j(15031) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(137) ! ALD4OH -> NII products
      j(15032) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(144) ! C5nALDOH -> NI products
      j(15041) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(162) ! C6nALDOH -> NI products
      j(15051) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(170) ! C7nALDOH -> NI products
      j(15061) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(215) ! C8nALDOH -> NI products
      j(15071) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(185) ! intAldOH -> R + CHO
      j(15111) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(175) ! Glycidaldehyde -> oxyranyl radical + CHO
      j(16021) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(65) ! CH3COCH3 -> CH3CO + CH3
      j(21011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(66) ! CH3COCH3 -> CO + 2 CH3
      j(21012) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(67) ! CH3COCH2CH3 -> CH3CO + CH2CH3
      j(21021) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(230) ! C3H7COCH3 -> CH3CO + C3H7
      j(21031) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(232) ! C3H7COCH3 -> C3H7 + CO + CH3
      j(21032) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(233) ! C3H7COCH3 -> CH3C(OH)=CH2 + CH2=CH2
      j(21033) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(229) ! C2H5COC2H5 -> C2H5CO + C2H5
      j(21041) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(234) ! C4H9COCH3 -> CH3CH=CH2 + CH2=C(OH)CH3
      j(21051) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(235) ! C3H7COC2H5 -> C2H5CO + C3H7
      j(21061) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(236) ! C3H7COC2H5 -> C3H7CO + C2H5
      j(21062) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(237) ! C3H7COC2H5 -> C3H7 + CO + C2H5
      j(21063) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(238) ! C3H7COC2H5 -> C2H5C(OH)=CH2 + CH2=CH2
      j(21064) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(364) ! lKET5 -> products
      j(21071) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(239) ! 4-heptanone -> NI products
      j(21081) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(241) ! 4-octanone -> NI products
      j(21091) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(260) ! MIPK -> CH3CO + i-C3H7
      j(22011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(261) ! MIPK -> i-C3H7CO + CH3
      j(22012) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(262) ! MIPK -> i-C3H7 + CO + CH3
      j(22013) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(263) ! MIPK -> CH2=CHOH + CH3CH=CH2
      j(22014) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(264) ! MIBK -> CH3CO + i-C4H9
      j(22021) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(270) ! 5-Me-2-hexanone -> CH3CO + CH2CH2CH(CH3)2
      j(22031) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(268) ! 4-Me-2-hexanone -> CH3C(OH)=CH2 + 2-butene
      j(22041) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(269) ! 4-Me-2-hexanone -> CH3C(OH)=CH2 + 1-butene
      j(22042) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(274) ! CH3CH(CH3)COCH(CH3)2 -> i-C3H7CO + i-C3H7
      j(22051) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(248) ! 2-Me-4-heptanone -> NI products
      j(22061) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(251) ! 3-Me-4-heptanone -> NI products
      j(22071) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(253) ! 2,2-Me-3-hexanone -> NI products
      j(22081) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(255) ! DIBK -> NI products
      j(22091) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(257) ! di-sec-butyl ketone -> NI products
      j(22101) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(259) ! di-t-butyl ketone -> NI products
      j(22111) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(244) ! n-C3H7COCH(CH3)2 -> NI products
      j(22121) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(245) ! n-C3H7COCH(CH3)2 -> i-C3H7COCH3 + C2H4
      j(22131) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(246) ! CH3COCH2C(CH3)3 -> NI products
      j(22141) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(258) ! di-sec-butyl ketone -> sec-C4H9COCH2CH3 + C2H4
      j(22151) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(59) ! CH3COCH=CH2 -> CH3 + C2H3CO
      j(23011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(295) ! CH3CH2COCH=CH2 -> C2H5 + C2H3CO
      j(23021) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(276) ! c-C3H4O -> C2H4 + CO
      j(24011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(277) ! c-C3H4O -> further products
      j(24012) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(279) ! c-C4H6O -> C2H4 + CH2=C=O
      j(24021) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(281) ! c-C4H6O -> c-C3H6 + CO
      j(24022) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(283) ! c-C5H8O -> 2 C2H4 + CO
      j(24031) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(284) ! c-C5H8O -> c-C4H8 + CO
      j(24032) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(285) ! c-C5H8O -> CH2=CHCH2CH2CHO
      j(24033) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(287) ! c-C6H10O -> 5-hexenal
      j(24041) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(289) ! c-C6H10O -> 1-pentene + CO
      j(24042) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(291) ! c-C7H12O -> 6-heptenal
      j(24051) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(293) ! c-C7H12O -> 1-hexene + CO
      j(24052) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(68) ! CH2(OH)COCH3 -> CH3CO + CH2(OH)
      j(25011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(299) ! CH3COC2H4OH -> CH3 + COCH2CH2OH
      j(25021) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(301) ! CH3COCH(OH)CH3 -> CH3CO + CH3CHOH
      j(25031) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(303) ! CH3COC(CH3)2OH -> CH3 + (CH3)2C(OH)CO
      j(25041) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(79) ! CH3COCOOH -> CH3CHO + CO2
      j(25071) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(304) ! CH3COCH2C(CH3)2OH -> CH3COCH2 + CH3C(OH)CH3
      j(25081) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(310) ! CH2=C=O -> CO2 + CO + H2
      j(26011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(311) ! CH3CH=C=O -> C2H4 + CO
      j(26012) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(70) ! CHOCHO -> 2 HO2 + 2 CO
      j(31011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(72) ! CHOCHO -> CH2O + CO
      j(31012) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(73) ! Ald (mult)
      j(31013) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(75) ! CH3COCHO -> CH3CO + HCO
      j(31021) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(76) ! CH3COCOCH3 -> Products
      j(31031) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(313) ! CHOCH=CHCHO -> 3H-furan-2-one
      j(33011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(315) ! CH3COCH=CHCHO -> 5Me-3H-2-furanone
      j(33021) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(316) ! CH3COCH=CHCHO -> CH3 + CHOCH=CHCO
      j(33022) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(317) ! CH3COCH=CHCHO -> CH3COCH=CH2 + CO
      j(33023) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(321) ! CH3COCH=CHCOCH3 -> CH3CO + CH=CHCOCH3
      j(33031) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(320) ! CHOCH=CHCH=CHCHO -> diformyl cyclobutene
      j(33041) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(322) ! pinonaldehyde -> R + CO + HO2
      j(34011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(323) ! caronaldehyde -> R + CO + HO2
      j(34021) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(34) ! CH3ONO2 -> CH3O + NO2
      j(41011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(37) ! C2H5ONO2 -> C2H5O + NO2
      j(41021) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(38) ! n-C3H7ONO2 -> C3H7O + NO2
      j(41031) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(39) ! 1-C4H9ONO2 -> 1-C4H9O + NO2
      j(41041) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(324) ! n-C5H11ONO2 -> n-C5H11O + NO2
      j(41051) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(41) ! CH3CHONO2CH3 -> CH3CHOCH3 + NO2
      j(41111) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(40) ! 2-C4H9ONO2 -> 2-C4H9O + NO2
      j(41121) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(325) ! 2-C5H11ONO2 -> 2-C5H11O + NO2
      j(41131) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(326) ! 3-C5H11ONO2 -> 3-C5H11O + NO2
      j(41141) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(327) ! C5H11ONO2 -> C5H11O + NO2
      j(41151) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(329) ! i-C4H9ONO2 -> i-C4H9O + NO2
      j(42011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(44) ! C(CH3)3(ONO2) -> C(CH3)3(O.) + NO2
      j(42021) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(330) ! i-C5H11ONO2 -> i-C5H11O + NO2
      j(42031) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(328) ! c-C5H11ONO2 -> c-C5H11O + NO2
      j(44011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(42) ! CH2(OH)CH2(ONO2) -> CH2(OH)CH2(O.) + NO2
      j(45011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(331) ! C1(OH)NO3 -> C1(OH)O + NO2
      j(45021) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(332) ! R(OH)NO3 -> R(OH)O + NO2
      j(45031) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(333) ! iR(OH)NO3 -> iR(OH)O + NO2
      j(45041) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(334) ! tR(OH)NO3 -> tR(OH)O + NO2
      j(45051) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(360) ! CH3NO2 -> CH3 + NO2
      j(46011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(361) ! RNO2 -> alkene + HONO
      j(46012) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(362) ! C2H5NO2 -> C2H5 + NO2
      j(46021) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(363) ! C2H5NO2 -> C2H4 + HONO
      j(46022) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(35) ! CH3(OONO2) -> CH3(OO) + NO2
      j(47011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(46) ! CH3CO(OONO2) -> CH3CO(OO) + NO2
      j(48011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(47) ! CH3CO(OONO2) -> CH3CO(O) + NO3
      j(48012) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(48) ! CH3CH2CO(OONO2) -> CH3CH2CO(OO) + NO2
      j(48021) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(49) ! CH3CH2CO(OONO2) -> CH3CH2CO(O) + NO3
      j(48022) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(335) ! PAN -> RCO(OO) + NO2
      j(48031) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(336) ! PAN -> RCO(O) + NO3
      j(48032) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(337) ! CH3CH(NO3)CH2NO3 -> CH3CH(NO3)CH2O + NO2
      j(51011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(339) ! C2H5CH(NO3)CH2NO3 -> C2H5CH(NO3)CH2O + NO2
      j(51021) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(341) ! CH3CH(NO3)CH(NO3)CH3 -> RO. + NO2
      j(51031) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(342) ! CH2(NO3)CH=CHCH2NO3 -> RO. + NO2
      j(53011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(343) ! CH2=CHCH(NO3)CH2NO3 -> CH2=CHCH(NO3)CH2O + NO2
      j(53021) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(32) ! CH3OOH -> CH3O + OH
      j(61011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(348) ! (CH3)3COOH -> (CH3)3CO + OH
      j(62011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(33) ! HOCH2OOH -> HOCH2O. + OH
      j(65011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(78) ! CH3CO(OOH) -> Products
      j(66011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(354) ! CH2OO -> HCHO + O(3P)
      j(71011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(355) ! CH3CHOO -> CH3CHO + O(3P)
      j(71021) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(356) ! synCH3CHOO -> CH3CHO + O(3P)
      j(71022) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(357) ! antiCH3CHOO -> CH3CHO + O(3P)
      j(71023) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(358) ! C2H5CHOO -> C2H5CHO + O(3P)
      j(71031) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(359) ! (CH3)2COO -> CH3COCH3 + O(3P)
      j(72011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(43) ! CH3COCH2(ONO2) -> CH3COCH2(O.) + NO2
      j(81011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(349) ! C2H5COCH2NO3 -> C2H5COCH2O + NO2
      j(81021) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(350) ! CH3COCH(NO3)CH3 -> CH3COCH(O.)CH3 + NO2
      j(81031) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(352) ! CH3COCH2CH2CH(OOH)CH3 -> RO. + OH
      j(81111) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(353) ! oxohexyl-hydroperoxide -> RO. + OH
      j(81121) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(27) ! genCH3CHO(poly)
      j(82011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(30) ! genC2H5CHO(poly)
      j(82021) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(138) ! C4nALDpoly
      j(82031) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(147) ! C5nALDpoly
      j(82041) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(165) ! C6nALDpoly
      j(82051) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(173) ! C7nALDpoly
      j(82061) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(218) ! nALDpoly(C>7)
      j(82071) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(226) ! genluALD(poly)
      j(82181) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(201) ! genbMeuAld(poly)
      j(82121) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(193) ! genaMeuAld(poly)
      j(82131) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(209) ! genluuALD(poly)
      j(82211) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(31) ! genC2H5CHO(OHpoly)
      j(83011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(139) ! C4nALDOHpoly
      j(83021) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(148) ! C5nALDOHpoly
      j(83031) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(166) ! C6nALDOHpoly
      j(83041) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(174) ! C7nALDOHpoly
      j(83051) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(219) ! nALDOHpoly(C>7)
      j(83061) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(227) ! genluALD(OHpoly)
      j(83111) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(202) ! genbMeuAldOH(poly)
      j(83121) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(194) ! genaMeuAldOH(poly)
      j(83131) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(210) ! genluuALD(OHpoly)
      j(83211) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(278) ! genC3cKet
      j(84011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(282) ! genC4cKet
      j(84021) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(286) ! genC5cKet
      j(84031) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(290) ! genC6cKet
      j(84041) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(294) ! genC7cKet
      j(84051) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(298) ! genuKet(poly)
      j(85011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(312) ! genKete(poly)
      j(86011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(74) ! Ald (poly)
      j(88011) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(314) ! uDICARaa(poly)
      j(88111) = seval(szabin,theta,tmp,tmp2,b,c,d)
    CASE(319) ! uDICARak(poly)
      j(88121) = seval(szabin,theta,tmp,tmp2,b,c,d)
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
