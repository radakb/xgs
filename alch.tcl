# alch.tcl
#
# This file exposes XGS hooks for performing expanded alchemical ensemble
# simulations in NAMD.
#
source [file join [file dirname [info script]] "xgs/xgs.alch.tcl"]

namespace eval ::alch {
    namespace import ::xgs::*
    rename xgsParamLadder alchLambdaLadder
    rename xgsStateWeights alchLambdaWeights
    rename xgsGibbsMethod alchGibbsMethod
    rename xgsRestartFile alchRestartFile

    namespace export alchRun alchCalibrate\
            alchLambdaLadder alchLambdaWeights alchGibbsMethod\
            alchRestartFile alchOptimalLambdaCount alchOptimalLambdaLadder
}

namespace import ::alch::*

