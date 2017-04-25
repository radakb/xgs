# st.tcl
#
# This file exposes XGS hooks for performing simulated tempering in NAMD.
#
# Simulated tempering (ST) can be loosely thought of as a single replica 
# version of temperature replica exchange (i.e. parallel tempering). The main 
# variation is that a set of fixed temperature dependent weights must be 
# specified at the beginning of the simulation and the choice of these will 
# greatly affect the sampling efficiency.
#
source [file join [file dirname [info script]] "xgs/xgs.st.tcl"]
source [file join [file dirname [info script]] "xgs/xgs.util.tcl"]

namespace eval ::st {
    namespace import ::xgs::*
    namespace import ::xgs::util::geometricTemperatures \
            ::xgs::util::optimalTempLadder

    rename xgsRun stRun
    rename xgsParamLadder stTempLadder
    rename xgsStateWeights stTempWeights
    rename xgsGibbsMethod stGibbsMethod
    rename xgsRestartFile stRestartFile
    rename xgsRestartFreq stRestartFreq

    namespace export stRun stCalibrate \
            stTempLadder stTempWeights stGibbsMethod \
            stRestartFile stRestartFreq \
            geometricTemperatures optimalTempLadder
}

namespace import ::st::*

