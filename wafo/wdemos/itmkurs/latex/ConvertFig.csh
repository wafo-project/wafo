#! /bin/tcsh -f

#
# Convert fig-files to LaTeX+eps-files
#

# Chapter Introduction
fig2pstex fig/FigRFCdef_intro
#fig2pstex fig/FigIntroOverview

# Chapter Cycles
#fig2pstex fig/FigTP_IntCross
#fig2pstex fig/FigDefRFC
#fig2pstex fig/FigRFC_minMax_def
fig2pstex fig/FigTP_Matrix

# Chapter Rainflow Cycles

#fig2pstex fig/FigProofRFM

# Chapter Asymetric Rainflow Cycles

#fig2pstex fig/FigDefARFC
#fig2pstex fig/FigTP_ARFM
#fig2pstex fig/FigHanging
#fig2pstex fig/FigStanding
#fig2pstex fig/FigEvents


# Chapter Examples
#fig2pstex fig/FigFlowChart1
fig2pstex fig/FigFlowChart2

# Chapter Decomposition


# Chapter Exact Distribution

#fig2pstex fig/FigExactDistr_IntCross
#fig2pstex fig/FigExactDistr_IntCrossR
