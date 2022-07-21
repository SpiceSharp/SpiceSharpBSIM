using SpiceSharp;

namespace SpiceSharpBSIM.Components.Semiconductors.BSIM
{

    public static class ACM
    {
        public static void SourceDrainResistances(int ACM, double LD, double LDIF, double HDIF, double WMLT, double w, double XW, double RSH,
            double RD, double RDC, GivenParameter<double> drainSquares, double RS, double RSC, GivenParameter<double> sourceSquares,
            out double drainConductance, out double sourceConductance)
        {
            switch (ACM)
            {
                case 1:
                case 11:
                    drainConductance = (LD + LDIF) / (w * WMLT + XW) * RD + RSH * drainSquares + RDC;
                    sourceConductance = (LD + LDIF) / (w * WMLT + XW) * RS + RSH * sourceSquares + RSC;

                    break;

                case 2:
                case 12:
                case 3:
                case 13:
                    if (drainSquares.Given)
                        drainConductance = (LD + LDIF) / (w * WMLT + XW) * RD + RSH * drainSquares + RDC;
                    else
                        drainConductance = ((LD + LDIF) * RD + (HDIF * WMLT) * RSH) / (w * WMLT + XW) + RDC;
                    if (sourceSquares.Given)
                        sourceConductance = (LD + LDIF) / (w * WMLT + XW) * RS + RSH * sourceSquares + RSC;
                    else
                        sourceConductance = ((LD + LDIF) * RS + (HDIF * WMLT) * RSH) / (w * WMLT + XW) + RSC;

                    break;

                default:
                    drainConductance = 0;
                    sourceConductance = 0;
                    break;
            }
        }

        /* Area Calculation Method (ACM) for MOS models */
        public static void SaturationCurrents(int ACM, int CALCACM, int GEO, double HDIF, double WMLT, double w,
            double XW, double jctTempSatCurDensity, double jctSidewallTempSatCurDensity,
            GivenParameter<double> drainArea, GivenParameter<double> drainPerimeter, GivenParameter<double> sourceArea, GivenParameter<double> sourcePerimeter,
            out double DrainSatCurrent, out double SourceSatCurrent)
        {
            switch (ACM)
            {
                case 1:
                case 11:
                    drainArea = (w * WMLT + XW) * WMLT;
                    drainPerimeter = (w * WMLT + XW);
                    DrainSatCurrent = drainArea * jctTempSatCurDensity + drainPerimeter * jctSidewallTempSatCurDensity;
                    if (DrainSatCurrent <= 0.0) DrainSatCurrent = 1.0e-14;

                    sourceArea = (w * WMLT + XW) * WMLT;
                    sourcePerimeter = (w * WMLT + XW);
                    SourceSatCurrent = sourceArea * jctTempSatCurDensity + sourcePerimeter * jctSidewallTempSatCurDensity;
                    if (SourceSatCurrent <= 0.0) SourceSatCurrent = 1.0e-14;

                    break;

                case 2:
                case 12:
                    if ((ACM == 2) || ((ACM == 12) && (CALCACM == 1)))
                    {
                        if (!drainArea.Given)
                            drainArea = 2.0 * (HDIF * WMLT) * (w * WMLT + XW);
                        else
                            drainArea = drainArea * WMLT * WMLT;
                        if (!drainPerimeter.Given)
                            drainPerimeter = 4.0 * (HDIF * WMLT) + 2.0 * (w * WMLT + XW);
                        else
                            drainPerimeter = drainPerimeter * WMLT;
                    }
                    DrainSatCurrent = drainArea * jctTempSatCurDensity + drainPerimeter * jctSidewallTempSatCurDensity;
                    if (DrainSatCurrent <= 0.0) DrainSatCurrent = 1.0e-14;

                    if ((ACM == 2) || ((ACM == 12) && (CALCACM == 1)))
                    {
                        if (!sourceArea.Given)
                            sourceArea = 2.0 * (HDIF * WMLT) * (w * WMLT + XW);
                        else
                            sourceArea = sourceArea * WMLT * WMLT;
                        if (!sourcePerimeter.Given)
                            sourcePerimeter = 4.0 * (HDIF * WMLT) + 2.0 * (w * WMLT + XW);
                        else
                            sourcePerimeter = sourcePerimeter * WMLT;
                    }
                    SourceSatCurrent = sourceArea * jctTempSatCurDensity + sourcePerimeter * jctSidewallTempSatCurDensity;
                    if (SourceSatCurrent <= 0.0) SourceSatCurrent = 1.0e-14;

                    break;

                case 3:
                case 13:
                    if (!drainArea.Given)
                        if ((GEO == 0) || (GEO == 2))
                            drainArea = 2.0 * (HDIF * WMLT) * (w * WMLT + XW);
                        else
                            drainArea = (HDIF * WMLT) * (w * WMLT + XW);
                    else
                        drainArea = drainArea * WMLT * WMLT;
                    if (!drainPerimeter.Given)
                        if ((GEO == 0) || (GEO == 2))
                            drainPerimeter = 4.0 * (HDIF * WMLT) + (w * WMLT + XW);
                        else
                            drainPerimeter = 2.0 * (HDIF * WMLT);
                    else
                        drainPerimeter = drainPerimeter * WMLT;
                    DrainSatCurrent = drainArea * jctTempSatCurDensity + drainPerimeter * jctSidewallTempSatCurDensity;
                    if (DrainSatCurrent <= 0.0) DrainSatCurrent = 1.0e-14;

                    if (!sourceArea.Given)
                        if ((GEO == 0) || (GEO == 1))
                            sourceArea = 2.0 * (HDIF * WMLT) * (w * WMLT + XW);
                        else
                            sourceArea = (HDIF * WMLT) * (w * WMLT + XW);
                    else
                        sourceArea = sourceArea * WMLT * WMLT;
                    if (!sourcePerimeter.Given)
                        if ((GEO == 0) || (GEO == 1))
                            sourcePerimeter = 4.0 * (HDIF * WMLT) + (w * WMLT + XW);
                        else
                            sourcePerimeter = 2.0 * (HDIF * WMLT);
                    else
                        sourcePerimeter = sourcePerimeter * WMLT;
                    SourceSatCurrent = sourceArea * jctTempSatCurDensity + sourcePerimeter * jctSidewallTempSatCurDensity;
                    if (SourceSatCurrent <= 0.0) SourceSatCurrent = 1.0e-14;

                    break;

                default:
                    DrainSatCurrent = 0;
                    SourceSatCurrent = 0;
                    break;
            }
        }

        public static void JunctionCapacitances(int ACM, int CALCACM, int GEO,
            double HDIF, double WMLT, double w, double XW,
            GivenParameter<double> drainArea, GivenParameter<double> drainPerimeter, GivenParameter<double> sourceArea, GivenParameter<double> sourcePerimeter,
            double CJ, double CJSW, double CJGATE,
            out double areaDrainBulkCapacitance, out double periDrainBulkCapacitance, out double gateDrainBulkCapacitance,
            out double areaSourceBulkCapacitance, out double periSourceBulkCapacitance, out double gateSourceBulkCapacitance)
        {
            switch (ACM)
            {
                case 1:
                    drainArea = (w * WMLT + XW) * WMLT;
                    drainPerimeter = (w * WMLT + XW);
                    areaDrainBulkCapacitance = drainArea * CJ;
                    periDrainBulkCapacitance = drainPerimeter * CJSW;
                    gateDrainBulkCapacitance = 0.0;

                    sourceArea = (w * WMLT + XW) * WMLT;
                    sourcePerimeter = (w * WMLT + XW);
                    areaSourceBulkCapacitance = sourceArea * CJ;
                    periSourceBulkCapacitance = sourcePerimeter * CJSW;
                    gateSourceBulkCapacitance = 0.0;

                    break;

                case 2:
                    if (!drainArea.Given)
                        drainArea = 2.0 * (HDIF * WMLT) * (w * WMLT + XW);
                    else
                        drainArea = drainArea * WMLT * WMLT;
                    if (!drainPerimeter.Given)
                        drainPerimeter = 4.0 * (HDIF * WMLT) + 2.0 * (w * WMLT + XW);
                    else
                        drainPerimeter = drainPerimeter * WMLT;
                    areaDrainBulkCapacitance = drainArea * CJ;
                    if (drainPerimeter > (w * WMLT + XW))
                    {
                        periDrainBulkCapacitance = (drainPerimeter - (w * WMLT + XW)) * CJSW;
                        gateDrainBulkCapacitance = (w * WMLT + XW) * CJGATE;
                    }
                    else
                    {
                        periDrainBulkCapacitance = drainPerimeter * CJGATE;
                        gateDrainBulkCapacitance = 0.0;
                    }

                    if (!sourceArea.Given)
                        sourceArea = 2.0 * (HDIF * WMLT) * (w * WMLT + XW);
                    else
                        sourceArea = sourceArea * WMLT * WMLT;
                    if (!sourcePerimeter.Given)
                        sourcePerimeter = 4.0 * (HDIF * WMLT) + 2.0 * (w * WMLT + XW);
                    else
                        sourcePerimeter = sourcePerimeter * WMLT;
                    areaSourceBulkCapacitance = sourceArea * CJ;
                    if (sourcePerimeter > (w * WMLT + XW))
                    {
                        periSourceBulkCapacitance = (sourcePerimeter - (w * WMLT + XW)) * CJSW;
                        gateSourceBulkCapacitance = (w * WMLT + XW) * CJGATE;
                    }
                    else
                    {
                        periSourceBulkCapacitance = sourcePerimeter * CJGATE;
                        gateSourceBulkCapacitance = 0.0;
                    }

                    break;

                case 3:
                    if (!drainArea.Given)
                        if ((GEO == 0) || (GEO == 2))
                            drainArea = 2.0 * (HDIF * WMLT) * (w * WMLT + XW);
                        else
                            drainArea = (HDIF * WMLT) * (w * WMLT + XW);
                    else
                        drainArea = drainArea * WMLT * WMLT;
                    if (!drainPerimeter.Given)
                        if ((GEO == 0) || (GEO == 2))
                            drainPerimeter = 4.0 * (HDIF * WMLT) + (w * WMLT + XW);
                        else
                            drainPerimeter = 2.0 * (HDIF * WMLT);
                    else
                        drainPerimeter = drainPerimeter * WMLT;
                    areaDrainBulkCapacitance = drainArea * CJ;
                    periDrainBulkCapacitance = drainPerimeter * CJSW;
                    gateDrainBulkCapacitance = (w * WMLT + XW) * CJGATE;

                    if (!sourceArea.Given)
                        if ((GEO == 0) || (GEO == 1))
                            sourceArea = 2.0 * (HDIF * WMLT) * (w * WMLT + XW);
                        else
                            sourceArea = (HDIF * WMLT) * (w * WMLT + XW);
                    else
                        sourceArea = sourceArea * WMLT * WMLT;
                    if (!sourcePerimeter.Given)
                        if ((GEO == 0) || (GEO == 1))
                            sourcePerimeter = 4.0 * (HDIF * WMLT) + (w * WMLT + XW);
                        else
                            sourcePerimeter = 2.0 * (HDIF * WMLT);
                    else
                        sourcePerimeter = sourcePerimeter * WMLT;
                    areaSourceBulkCapacitance = sourceArea * CJ;
                    periSourceBulkCapacitance = sourcePerimeter * CJSW;
                    gateSourceBulkCapacitance = (w * WMLT + XW) * CJGATE;

                    break;

                case 11:
                    drainArea = (w * WMLT + XW) * WMLT;
                    drainPerimeter = (w * WMLT + XW);
                    areaDrainBulkCapacitance = drainArea * CJ;
                    periDrainBulkCapacitance = drainPerimeter * CJSW;
                    gateDrainBulkCapacitance = 0.0;

                    sourceArea = (w * WMLT + XW) * WMLT;
                    sourcePerimeter = (w * WMLT + XW);
                    areaSourceBulkCapacitance = sourceArea * CJ;
                    periSourceBulkCapacitance = sourcePerimeter * CJSW;
                    gateSourceBulkCapacitance = 0.0;

                    break;

                case 12:
                    if (CALCACM == 1)
                    {
                        if (!drainArea.Given)
                            drainArea = 2.0 * (HDIF * WMLT) * (w * WMLT + XW);
                        else
                            drainArea = drainArea * WMLT * WMLT;
                        if (!drainPerimeter.Given)
                            drainPerimeter = 4.0 * (HDIF * WMLT) + 2.0 * (w * WMLT + XW);
                        else
                            drainPerimeter = drainPerimeter * WMLT;
                    }
                    areaDrainBulkCapacitance = drainArea * CJ;
                    if (drainPerimeter > (w * WMLT + XW))
                    {
                        periDrainBulkCapacitance = (drainPerimeter - (w * WMLT + XW)) * CJSW;
                        gateDrainBulkCapacitance = (w * WMLT + XW) * CJGATE;
                    }
                    else
                    {
                        periDrainBulkCapacitance = 0.0;
                        gateDrainBulkCapacitance = drainPerimeter * CJGATE;
                    }

                    if (CALCACM == 1)
                    {
                        if (!sourceArea.Given)
                            sourceArea = 2.0 * (HDIF * WMLT) * (w * WMLT + XW);
                        else
                            sourceArea = sourceArea * WMLT * WMLT;
                        if (!sourcePerimeter.Given)
                            sourcePerimeter = 4.0 * (HDIF * WMLT) + 2.0 * (w * WMLT + XW);
                        else
                            sourcePerimeter = sourcePerimeter * WMLT;
                    }
                    areaSourceBulkCapacitance = sourceArea * CJ;
                    if (sourcePerimeter > (w * WMLT + XW))
                    {
                        periSourceBulkCapacitance = (sourcePerimeter - (w * WMLT + XW)) * CJSW;
                        gateSourceBulkCapacitance = (w * WMLT + XW) * CJGATE;
                    }
                    else
                    {
                        periSourceBulkCapacitance = 0.0;
                        gateSourceBulkCapacitance = sourcePerimeter * CJGATE;
                    }

                    break;

                case 13:
                    drainArea = drainArea * WMLT * WMLT;
                    drainPerimeter = drainPerimeter * WMLT;
                    areaDrainBulkCapacitance = drainArea * CJ;
                    if (drainPerimeter > (w * WMLT + XW))
                    {
                        periDrainBulkCapacitance = (drainPerimeter - (w * WMLT + XW)) * CJSW;
                        gateDrainBulkCapacitance = (w * WMLT + XW) * CJGATE;
                    }
                    else
                    {
                        periDrainBulkCapacitance = 0.0;
                        gateDrainBulkCapacitance = drainPerimeter * CJGATE;
                    }

                    sourceArea = sourceArea * WMLT * WMLT;
                    sourcePerimeter = sourcePerimeter * WMLT;
                    areaSourceBulkCapacitance = sourceArea * CJ;
                    if (sourcePerimeter > (w * WMLT + XW))
                    {
                        periSourceBulkCapacitance = (sourcePerimeter - (w * WMLT + XW)) * CJSW;
                        gateSourceBulkCapacitance = (w * WMLT + XW) * CJGATE;
                    }
                    else
                    {
                        periSourceBulkCapacitance = 0.0;
                        gateSourceBulkCapacitance = sourcePerimeter * CJGATE;
                    }

                    break;

                default:
                    areaDrainBulkCapacitance = 0;
                    periDrainBulkCapacitance = 0;
                    gateDrainBulkCapacitance = 0;
                    areaSourceBulkCapacitance = 0;
                    periSourceBulkCapacitance = 0;
                    gateSourceBulkCapacitance = 0;
                    break;
            }
        }
    }
}
