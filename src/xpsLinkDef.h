#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

//StatUtils.h
#pragma link C++ class TUnivariateTest+;
#pragma link C++ class TStudentTest+;
#pragma link C++ class TStat+;
#pragma link C++ class TMLMath+;

//XPSUtils.h
#pragma link C++ enum EPlotErrors;
#pragma link C++ class XBitSet+;
#pragma link C++ class XPlot+;

//XPSBase.h
#pragma link C++ class XIdxString+;
#pragma link C++ class XLdxString+;
#pragma link C++ class XFolder+;
#pragma link C++ class XSetting+;
#pragma link C++ class XTreeInfo+;
#pragma link C++ class XTreeHeader+;
#pragma link C++ class XTreeSet+;
#pragma link C++ class XAlgorithm+;
#pragma link C++ class XManager+;

//XPSSchemes.h
#pragma link C++ class XSchemeManager+;
#pragma link C++ class XLayoutTreeInfo+;
#pragma link C++ class XProbeTreeInfo+;
#pragma link C++ class XSchemeTreeInfo+;
#pragma link C++ class XUnitTreeInfo+;
#pragma link C++ class XGenomeTreeInfo+;
#pragma link C++ class XDNAChip+;
#pragma link C++ class XMicroArray+;
#pragma link C++ class XOligoArray+;
#pragma link C++ class XGeneChip+;
#pragma link C++ class XSNPChip+;
#pragma link C++ class XGenomeChip+;
#pragma link C++ class XExonChip+;
#pragma link C++ class XPosition+;
#pragma link C++ class XLayout+;
#pragma link C++ class XScheme+;
#pragma link C++ class XGCScheme+;
#pragma link C++ class XExonScheme+;
#pragma link C++ class XProbe+;
#pragma link C++ class XGCProbe+;
#pragma link C++ class XUnitID+;
#pragma link C++ class XUnit+;
#pragma link C++ class XGCUnit+;
#pragma link C++ class XExonUnit+;
#pragma link C++ class XAnnotation+;
#pragma link C++ class XTransAnnotation+;
#pragma link C++ class XGenomeAnnotation+;
#pragma link C++ class XExonAnnotation+;
#pragma link C++ class XProbesetAnnotation+;

//XPSDataTypes.h
#pragma link C++ class XDataTypeInfo+;
#pragma link C++ class XDataTypeList+;
#pragma link C++ class XDatabaseInfo+;
#pragma link C++ class XProjectInfo+;
#pragma link C++ class XAuthorInfo+;
#pragma link C++ class XLoginInfo+;
#pragma link C++ class XDatasetInfo+;
#pragma link C++ class XSourceInfo+;
#pragma link C++ class XArrayInfo+;
#pragma link C++ class XHybInfo+;
#pragma link C++ class XHybridizationList+;
#pragma link C++ class XSampleInfo+;
#pragma link C++ class XCellLineInfo+;
#pragma link C++ class XPrimaryCellInfo+;
#pragma link C++ class XTissueInfo+;
#pragma link C++ class XBiopsyInfo+;
#pragma link C++ class XTreatmentInfo+;
#pragma link C++ class XTreatmentList+;

//XProjectHandler.h
#pragma link C++ class XHandler+;
#pragma link C++ class XProjectHandler+;

//XPSData.h
#pragma link C++ class XDataManager+;
#pragma link C++ class XDataTreeInfo+;
#pragma link C++ class XMaskTreeInfo+;
#pragma link C++ class XDataSetting+;
#pragma link C++ class XHybridization+;
#pragma link C++ class XGeneChipHyb+;
#pragma link C++ class XGeneChipMetrics+;
#pragma link C++ class XGeneChipPivot+;
#pragma link C++ class XSNPChipHyb+;
#pragma link C++ class XSNPChipMetrics+;
#pragma link C++ class XSNPChipPivot+;
#pragma link C++ class XGenomeChipHyb+;
#pragma link C++ class XGenomeChipMetrics+;
#pragma link C++ class XGenomeChipPivot+;
#pragma link C++ class XExonChipHyb+;
#pragma link C++ class XExonChipMetrics+;
#pragma link C++ class XExonChipPivot+;
#pragma link C++ class XGenePixHyb+;
#pragma link C++ class XCellOM+;
#pragma link C++ class XMask+;
#pragma link C++ class XCellMask+;
#pragma link C++ class XUnitMask+;
#pragma link C++ class XCell+;
#pragma link C++ class XGCCell+;
#pragma link C++ class XBgCell+;
#pragma link C++ class XRankCell+;
#pragma link C++ class XResidual+;
#pragma link C++ class XBorder+;
#pragma link C++ class XExpression+;
#pragma link C++ class XGCExpression+;
#pragma link C++ class XQCExpression+;
#pragma link C++ class XSpliceExpression+;
#pragma link C++ class XCall+;
#pragma link C++ class XPCall+;
#pragma link C++ class XGPCell+;
#pragma link C++ class XFeature635+;
#pragma link C++ class XFeature532+;
#pragma link C++ class XBg532+;
#pragma link C++ class XGPRatio+;
#pragma link C++ class XAnalysisPlot+;
#pragma link C++ class XGeneChipPlot+;
#pragma link C++ class XGenePixPlot+;

//XPSProcessing.h
#pragma link C++ class XProcessManager+;
#pragma link C++ class XExpressionTreeInfo+;
#pragma link C++ class XSelectionTreeInfo+;
#pragma link C++ class XProcesSetting+;
#pragma link C++ class XProcesSet+;

//XPSPreProcessing.h
#pragma link C++ class XPreProcessManager+;
#pragma link C++ class XBorderTreeInfo+;
#pragma link C++ class XCallTreeInfo+;
#pragma link C++ class XQualityTreeInfo+;
#pragma link C++ class XResidualTreeInfo+;
#pragma link C++ class XPreProcesSetting+;
#pragma link C++ class XPreProcesSet+;
#pragma link C++ class XGCProcesSet+;
#pragma link C++ class XGenomeProcesSet+;
#pragma link C++ class XExonProcesSet+;

//XPSSelector.h
#pragma link C++ class XSelector+;
#pragma link C++ class XRankSelector+;
#pragma link C++ class XProbeSelector+;
#pragma link C++ class XUnitSelector+;
#pragma link C++ class XUserSelector+;

//XPSSplicer.h
//#pragma link C++ class XSplicer+;
//#pragma link C++ class XSpliceIndex+;
//#pragma link C++ class XFIRMA+;

//XPSHybridizer.h
#pragma link C++ class XHybridizer+;
#pragma link C++ class XBackgrounder+;
#pragma link C++ class XCallDetector+;
#pragma link C++ class XExpressor+;
#pragma link C++ class XMultichipExpressor+;
#pragma link C++ class XSpliceExpressor+;
#pragma link C++ class XQualifier+;
#pragma link C++ class XSectorBackground+;
#pragma link C++ class XWeightedBackground+;
#pragma link C++ class XRMABackground+;
#pragma link C++ class XGCBackground+;
#pragma link C++ class XMeanDifferenceCall+;
#pragma link C++ class XDetectionCall+;
#pragma link C++ class XMAS4Call+;
#pragma link C++ class XDABGCall+;
#pragma link C++ class XINICall+;
#pragma link C++ class XArithmeticMean+;
#pragma link C++ class XGeometricMean+;
#pragma link C++ class XWeightedMean+;
#pragma link C++ class XGCCorrectedMean+;
#pragma link C++ class XWeightedDiff+;
#pragma link C++ class XAvgDif+;
#pragma link C++ class XTukeyBiweight+;
#pragma link C++ class XMedianPolish+;
#pragma link C++ class XFARMS+;
#pragma link C++ class XDFW+;
#pragma link C++ class XFIRMA+;
#pragma link C++ class XRMAQualifier+;
#pragma link C++ class XPLMQualifier+;

//XPSNormation.h
#pragma link C++ class XNormationManager+;
#pragma link C++ class XNormationSetting+;
#pragma link C++ class XNormedSet+;
#pragma link C++ class XNormedGCSet+;
#pragma link C++ class XNormedGenomeSet+;
#pragma link C++ class XNormedExonSet+;

//XPSNormalizer.h
#pragma link C++ class XNormalizer+;
#pragma link C++ class XMeanNormalizer+;
#pragma link C++ class XMedianNormalizer+;
#pragma link C++ class XKernelNormalizer+;
#pragma link C++ class XLowessNormalizer+;
#pragma link C++ class XSuperNormalizer+;
#pragma link C++ class XQuantileNormalizer+;

//XPSAnalysis.h
#pragma link C++ class XAnalysisManager+;
#pragma link C++ class XAnalySetting+;
#pragma link C++ class XPreFilterHeader+;
#pragma link C++ class XUniFilterHeader+;
#pragma link C++ class XMultiFilterHeader+;
#pragma link C++ class XUniTestHeader+;
#pragma link C++ class XMultiTestHeader+;
#pragma link C++ class XClusterHeader+;
#pragma link C++ class XRegressionHeader+;
#pragma link C++ class XAnalySet+;
#pragma link C++ class XPreFilterSet+;
#pragma link C++ class XUnivarSet+;
#pragma link C++ class XMultivarSet+;
#pragma link C++ class XClusterSet+;
#pragma link C++ class XRegressionSet+;
#pragma link C++ class XScore+;
#pragma link C++ class XGrpMn+;
#pragma link C++ class XChance+;
#pragma link C++ class XAdjust+;

//XPSAnalyzer.h
#pragma link C++ class XAnalyser+;
#pragma link C++ class XUniTester+;
#pragma link C++ class XMultiTester+;
#pragma link C++ class XClusterizer+;
#pragma link C++ class XRegressor+;

//XPSFilter.h
#pragma link C++ class XFilter+;
#pragma link C++ class XPreFilter+;
#pragma link C++ class XUniFilter+;
#pragma link C++ class XMultiFilter+;

#endif
