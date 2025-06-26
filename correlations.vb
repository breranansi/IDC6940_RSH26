'======================================
'== Molecular Weight Correlations    ==
'== Mw = f(tb,sg)                    ==
'======================================

Function hoss_Mw(tb As Double, sg As Double) As Double
'------------------------------------------------------------------------------
'-- Molecular Mass of Petroleum Fraction [lb/lb-mol]
'-- Hosseinifar & Shahverdi, 2021 (JPSE)
'-- Coded by Richard Henry, June 2025
'------------------------------------------------------------------------------
'-- tb = normal boiling point, [K]
'-- sg = specific gravity [1/wtr]
'------------------------------------------------------------------------------
    Const a00 As Double = 111.10258
    Const a01 As Double = -1.318569
    Const a02 As Double = -3.288424
    Const a03 As Double = -41.663306
    Const a04 As Double = -1.153088
    Const a05 As Double = -3.698071
    Const a06 As Double = -0.739074
    '--------------------------------------------------------------------------
    Mw = (3# + 2# * sg) / (3# - 2# * sg)
    hoss_Mw = a00 * (tb ^ a01) * (Mw ^ (a02 / 2#))
    hoss_Mw = hoss_Mw + a03 * (tb ^ a04) * (Mw ^ (a05 / 2#))
    hoss_Mw = hoss_Mw ^ a06
End Function

Function hari_Mw(tb As Double, sg As Double) As Double
'------------------------------------------------------------------------------
'-- Molecular Mass of Petroleum Fraction [lb/lb-mol]
'-- Hariu & Sage, 1969
'-- Coded by Richard Henry, June 2025
'------------------------------------------------------------------------------
'-- tb = normal boiling point, [K]
'-- sg = specific gravity [1/wtr]
'------------------------------------------------------------------------------
    Const a00 As Double = 0.6670202
    Const a01 As Double = 0.1552531
    Const a02 As Double = -0.005378496
    Const a03 As Double = 0.004583705
    Const a04 As Double = 0.00002500584
    Const a05 As Double = 0.000002698693
    Const a06 As Double = 0.0000387595
    Const a07 As Double = -0.00000001566228
    '--------------------------------------------------------------------------
    kw = ((1.8 * tb) ^ (1 / 3)) / sg
    Mw = a00 + a01 * kw + a02 * (kw ^ 2) + a03 * (tb * kw)
    Mw = Mw + a04 * (tb * (kw ^ 2)) + a05 * (tb ^ 2) + a06 * ((tb ^ 2) * kw)
    Mw = Mw + a07 * (tb ^ 2) * (kw ^ 2)
    hari_Mw = Mw
End Function

Function lina_Mw(tb As Double, sg As Double) As Double
'------------------------------------------------------------------------------
'-- Molecular Mass of Petroleum Fraction [lb/lb-mol]
'-- Correlation for Predicting the Molecular Weight of Brazilian Petroleum
'-- Residues and Cuts: Linan et al, 2011
'-- Coded by Richard Henry, June 2025
'------------------------------------------------------------------------------
'-- tb = normal boiling point, [K]
'-- sg = specific gravity [1/wtr]
'------------------------------------------------------------------------------
    Const a00 As Double = 284.75
    Const a01 As Double = 0.00322
    Const a02 As Double = -2.52
    Const a03 As Double = 0.083
    Const a04 As Double = 2.44
    '--------------------------------------------------------------------------
    lina_Mw = a00 * (Math.Exp(a01 * tb)) * (Math.Exp(a02 * sg))
    lina_Mw = lina_Mw * (tb ^ a03) * (sg ^ a04)
End Function

Function goos_Mw(tb As Double, sg As Double) As Double
'------------------------------------------------------------------------------
'-- Molecular Weight of a petroleum fraction [lb/lb-mol]                     --
'-- Prediction of Molecular Weight of Petroleum Fractions                    --
'-- Adriaan Goossens, Ind. Eng. Chem. Res. 1996, 35, 985-988                 --
'-- Coded by Richard Henry, March 2012                                       --
'------------------------------------------------------------------------------
'-- tb = normal boiling point, [degF] (91.4-1364.0)                          --
'-- sg = specific gravity [1/wtr] (0.63-1.08)                                --
'------------------------------------------------------------------------------
    Const c0 As Double = 0.01077
    Const c1 As Double = 1.52869
    Const c2 As Double = 0.06486
    Const c3 As Double = 1078#
    Dim c4 As Double
    Dim t As Double
    Dim d As Double
    Const c5 As Double = 1.002733049
    Const c6 As Double = -0.006433370483
    '--------------------------------------------------------------------------
    t = (tb + 459.68) / 1.8
    d = sg * c5 + c6
    c4 = c1 + c2 * Math.Log(t / (c3 - t))
    goos_Mw = c0 * (t ^ c4) / d
End Function

Function ria1_Mw(ByVal tb As Double, ByVal sg As Double) As Double
'-------------------------------------
'-- Molecular Weight [lb/mol]       --
'-- Later Riazi & Daubert First Form--
'-- Richard Henry October 2010      --
'-------------------------------------
'-- tb = normal boiling point [degF]--
'-- sg = specific gravity [1/water] --
'-------------------------------------
    Dim a As Double, b As Double, c As Double, d As Double, e As Double, f As Double
    Dim t As Double
    '-------------------------------------
    a = 581.96
    b = 0.97476
    c = 6.51274
    d = 0.000543076
    e = -9.53384
    f = 0.00111056
    '-------------------------------------
    t = tb + 460#
    ria1_Mw = a * (t ^ b) * (sg ^ c) * Math.Exp(d * t + e * sg + f * t * sg)
End Function

Function riaa_Mw(tb As Double, sg As Double) As Double
'-------------------------------------
'-- Molecular Weight [lb/mol]       --
'-- Early Riazi & Daubert           --
'-- Richard Henry October 2010      --
'-------------------------------------
'-- tb = normal boiling point [degF]--
'-- sg = specific gravity [1/water] --
'-------------------------------------
    Dim a As Double, b As Double, c As Double
    '-------------------------------------
    a = 0.000045673
    b = 2.1962
    c = -1.0164
    riaa_Mw = a * ((tb + 460#) ^ b) * (sg ^ c)
End Function

Function kesl_Mw(ByVal tb As Double, ByVal sg As Double) As Double
'-------------------------------------
'-- Molecular Weight [lb/mol]       --
'-- Kesler & Lee [1976]             --
'-- Richard Henry October 2010      --
'-------------------------------------
'-- tb = normal boiling point [degF]--
'-- sg = specific gravity [1/water] --
'-------------------------------------
    Dim a0 As Double, a1 As Double, a2 As Double, a3 As Double, a4 As Double, a5 As Double
    Dim a6 As Double, a7 As Double, a8 As Double, a9 As Double, aa As Double, ab As Double
    Dim t As Double
    '-------------------------------------
    a0 = -12272.6
    a1 = 9486.4
    a2 = 4.6523
    a3 = -3.3287
    a4 = -0.77084
    a5 = -0.02058
    a6 = 1.3437
    a7 = -720.79
    a8 = -0.80882
    a9 = 0.02226
    aa = 1.8828
    ab = -181.98
    '-------------------------------------
    t = tb + 460#
    kesl_Mw = a0 + a1 * sg + (a2 + a3 * sg) * t
    kesl_Mw = kesl_Mw + (1# + a4 * sg + a5 * sg ^ 2) * (a6 + a7 / t) * (10000000# / t)
    kesl_Mw = kesl_Mw + (1# + a8 * sg + a9 * sg ^ 2) * (aa + ab / t) * (1000000000000# / (t ^ 3))
End Function

Function api_Mw(tb As Double, sg As Double) As Double
'-------------------------------------
'-- Molecular Weight [lb/mol]       --
'-- American Petroleum Institute    --
'-- Whitson Page 79                 --
'-- Richard Henry March 2012        --
'-------------------------------------
'-- tb = normal boiling point [degF]--
'-- sg = specific gravity [1/water] --
'-------------------------------------
    Const a0 As Double = 204.38
    Const a1 As Double = 0.118
    Const a2 As Double = 1.88
    Const a3 As Double = 0.00218
    Const a4 As Double = -3.07
    Dim t As Double
    '-------------------------------------
    t = tb + 460#
    api_Mw = a0 * (t ^ a1) * (sg ^ a2) * Math.Exp(a3 * t + a4 * sg)
End Function

Function raob_Mw(tb As Double, sg As Double) As Double
'-------------------------------------
'-- Molecular Weight [lb/mol]       --
'-- Rao & Bardon (1985)             --
'-- Whitson Page 79                 --
'-- Richard Henry March 2012        --
'-------------------------------------
'-- tb = normal boiling point [degF]--
'-- sg = specific gravity [1/water] --
'-------------------------------------
    Const a0 As Double = 1.27
    Const a1 As Double = 0.071
    Const a2 As Double = 1.8
    Const a3 As Double = 22.31
    Const a4 As Double = 1.68
    Dim k As Double
    '--------------------------------------------------------------------------
    '-- Adjusted June 2025                                                   --
    '--------------------------------------------------------------------------
    k = (tb + 460#) ^ (1# / 3#)
    k = k / sg
    'rao_Mw = (a0 + a1 * k) * Math.Log(a2 * (Tb + 460#) / (a3 + a4 * k))
    raob_Mw = (a0 + a1 * k) * Math.Log((tb + 460#) / a2 / (a3 + a4 * k))
    raob_Mw = Math.Exp(raob_Mw)
End Function

Function wins_Mw(tb As Double, sg As Double) As Double
'-------------------------------------
'-- Molecular Weight [lb/mol]       --
'-- Winn(1957) Sim & Daubert (1980) --
'-- Richard Henry October 2010      --
'-------------------------------------
'-- tb = normal boiling point [degF]--
'-- sg = specific gravity [1/water] --
'-------------------------------------
    Dim a0 As Double, a1 As Double, a2 As Double
    '-------------------------------------
    a0 = 1.4350476 * 10 ^ -5
    a1 = 2.3776
    a2 = -0.9371
    '-------------------------------------
    wins_Mw = a0 * ((tb + 460#) ^ a1) * (sg ^ a2)
End Function

Function pede_Mw(tb As Double, sg As Double) As Double
'---------------------------------------
'-- Molecular weight [lb/mol]         --
'-- Pedersen, from Boiling Point      --
'-- Coded by Richard Henry, March 2012--
'---------------------------------------
'-- tb = Normal Boiling Point [degF]  --
'-- sg = Specific Gravity [1/air]     --
'---------------------------------------
    Const a As Double = 97.58
    Const b As Double = 0.3323
    Const c As Double = 0.4609
    '---------------------------------------
    pede_Mw = (tb + 459.68) / 1.8
    pede_Mw = (pede_Mw / (a * sg ^ c)) ^ (1# / b)
End Function

'======================================

function shossTb(sg as double) as double
'------------------------------------------------------------------------------
'-- Boiling Point [K] 
'------------------------------------------------------------------------------
'-- sg = Specific Gravity [1/air]     --
'------------------------------------------------------------------------------
const a00 as double = 50.810766
const a01 as double = 8.9376090
const a02 as double = 56.359421
const a03 as double = 0.851344
'------------------------------------------------------------------------------
tb=(3.0+2.0*sg)/(3.0-sg)
tb=a00*(tb ^ (a01/2.0))
tb=(tb+a02) ^ a03
shosstb=tb
end function

function mhossTb(mw as double) as double
'------------------------------------------------------------------------------
'-- Boiling Point [K] 
'------------------------------------------------------------------------------
'-- mw = Molecular Mass [lb/mol]
'------------------------------------------------------------------------------
const a00 as double = 0.440056
const a01 as double = -0.281715
const a02 as double = 0.288766
const a03 as double = -6.616655
'------------------------------------------------------------------------------
tb=a00*(mw ^ a01)
tb=(tb+a02) ^ a03
mhosstb=tb
end function
