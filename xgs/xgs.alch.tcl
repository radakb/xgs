# xgs.alch.tcl
#
# This file implements XGS protocols for running alchemical expanded ensemble 
# (ALCH) simulations. Rather than require NAMD users to handle the ::xgs 
# namespace directly, only those parts that are relevant to the user are 
# imported/exported in alch.tcl. 
#
package require Tcl 8.5

source [file join [file dirname [info script]] "xgs.core.tcl"]
source [file join [file dirname [info script]] "xgs.util.tcl"]

namespace eval ::xgs {
    variable TYPE ALCH
    namespace import ::xgs::util::linearLambdas \
            ::xgs::util::optimalLambdaCount ::xgs::util::optimalLambdaLadder
    namespace export alchCalibrate linearLambdas optimalLambdaCount \
            optimalLambdaLadder
}

# Here we assume that the second derivative of the free energy is a constant
# with respect to lambda. This is exactly true for a perturbation of a central
# potential with spherical symmetry (e.g. charging of a spherical ion, ala
# the Born equation). Two routes are available:
#
# 1) compute and average the numerical derivative of mean work between lambda
# values (this is a scaled average of the variances)
#
# 2) perform a linear regression of the form:
#   df(l) = 0.5*df2dl2*(1 - (l-1)**2)
# NB: This depends on the direction of the perturbation.
#
proc ::xgs::alchCalibrate {numsteps {numEquilSteps 0}} {
    lassign [xgsCalibrate $numsteps $numEquilSteps] \
            weightList varList energyMeans
    set lambdaList [getParams]
    set liList [lrange $lambdaList 0 end-1]
    set lip1List [lrange $lambdaList 1 end] 
    set lmin [lindex $lambdaList 0]
    set lmax [lindex $lambdaList end]
    # Compute the average second derivative on each interval.
    #
    set d2fdl2_mean 0.0
    foreach var $varList li $liList lip1 $lip1List { 
        set dl [expr {$lip1 - $li}]
        set d2fdl2_mean [expr {$d2fdl2_mean + $var/$dl**2}]
    }
    set d2fdl2_mean [expr {$d2fdl2_mean / [llength $varList]}]

    # Try a linear fit to du/dl instead.
    #
    set dudl [list]
    foreach meanList $energyMeans {
        set U1 [expr {lsum([lrange $meanList 0 2])}]
        set U2 [expr {lsum([lrange $meanList 3 end])}]
        set du [expr {($U1 - $U2) / ($::BOLTZMANN*[$::thermostatTempCmd])}]
        lappend dudl $du
    }
    lassign [linearFit $lambdaList $dudl] d2fdl2_fit1 U0 d2fdl2_err1 U0_err R1
    set d2fdl2_fit1 [expr {abs($d2fdl2_fit1)}]

    # Can also try to the fit free energy directly as linearized quadratic.
    # NB: The intercept is exactly 0.0 by construction - omit the reference
    #   free energy at (0, 0) and constrain the fit.
    #
    set x [list]
    if {[lindex $weightList end] < 0.0} { ;# for "charging" - dG < 0
        foreach lambda [lrange $lambdaList 1 end] {
            lappend x [expr {0.5*$lambda**2}]
        }
    } else { ;# for "decharging" - dG > 0
        foreach lambda [lrange $lambdaList 1 end] {
            lappend x [expr {0.5*(1 - (1 - $lambda)**2)}]
        }
    }
    set y [lrange $weightList 1 end]
    set y_err [lrepeat [llength $y] 1.0] ;# dummy values for weighted fit
    lassign [linearFit $x $y $y_err 0.0] d2fdl2_fit2 tmp d2fdl2_err2 tmp_err R2
    set d2fdl2_fit2 [expr {abs($d2fdl2_fit2)}]

    # Report ladders for each estimate above.
    #
    set N [optimalLambdaCount $d2fdl2_mean $lmin $lmax]
    set lladder [optimalLambdaLadder $d2fdl2_mean $lmin $lmax]
    xgsPrint "Properties from linear response:"
    xgsPrint "Using mean numerical derivative:"
    xgsPrint [format "|d2f/dl2| (kBT units): %11.4f" $d2fdl2_mean]
    xgsPrint "optimal lambda count: $N"
    xgsPrint "optimal lambda ladder: $lladder"
    set N [optimalLambdaCount $d2fdl2_fit1 $lmin $lmax]
    set lladder [optimalLambdaLadder $d2fdl2_fit1 $lmin $lmax]
    xgsPrint [format "Using linear fit to dU/dlambda (r^2 = %5.3f):" $R1]
    xgsPrint [format "|d2f/dl2| (kBT units): %11.4f +/- %11.4f" \
           $d2fdl2_fit1 $d2fdl2_err1]
    xgsPrint [format "intercept: %11.4f +/- %11.4f" $U0 $U0_err]
    xgsPrint "optimal lambda count: $N"
    xgsPrint "optimal lambda ladder: $lladder"
    set N [optimalLambdaCount $d2fdl2_fit2 $lmin $lmax]
    set lladder [optimalLambdaLadder $d2fdl2_fit2 $lmin $lmax]
    xgsPrint [format "Using linearized fit to weights (r^2 = %5.3f):" $R2]
    xgsPrint [format "|d2f/dl2| (kBT units): %11.4f +/- %11.4f" \
           $d2fdl2_fit2 $d2fdl2_err2]
    xgsPrint [format "(fixed) intercept: %11.4f +/- %11.4f" $tmp $tmp_err]
    xgsPrint "optimal lambda count: $N"
    xgsPrint "optimal lambda ladder: $lladder"
    return
}

# =============================================================================
# Gibbs Sampling Routines - necessary to run xgsRun
# =============================================================================
# ::xgs::parameterSetup
#
# Perform parameter setup specific to this method. This can also involve error
# checking that the user has not already set invalid parameters in the 
# configuration file.
#
proc ::xgs::parameterSetup {} {
    # The reduced potential routines are currently defined under the assumption
    # of strictly linear alchemical coupling.
    #
    if {[alchVdwShiftCoeff] > 0. || ![isset alchVdwShiftCoeff]} {
        xgsAbort "Use with soft-core potentials is not currently supported."
    }

    set Lambda [getParams [getStateIndex]]
    alchLambda $Lambda
    setStateIndex [getParamIndex $Lambda]
    return
}

# ::xgs::updateState
#
# Update the parameters after a successful state change. 
#
# Arguments:
# ----------
# oldIndex : int
#   The index of the old state
# newIndex : int
#   The index of the new (accepted) state
#
# Returns:
# --------
# None 
# 
proc ::xgs::updateState {oldIndex newIndex} {
    alchLambda [getParams $newIndex]
    return
}

# ::xgs:computePairPotential
#
# Return the reduced potential difference between two states.
#
# Arguments:
# ----------
# oldIndex : int
#   The current state index
# newIndex : int
#   The proposed state index
#
# Returns:
# --------
# du : float 
#   The reduced potential energy difference 
# 
proc ::xgs::computePairPotential {oldIndex newIndex} {
    set oldLambda [getParams $oldIndex]
    set newLambda [getParams $newIndex]
    set dweight [expr {[getWeights $newIndex] - [getWeights $oldIndex]}]
    set dU 0.0
    foreach li [getLambdas $oldLambda] lj [getLambdas $newLambda]\
            U [getAlchEnergies ::energyArray] {
        set dU [expr {$dU + ($lj - $li)*$U}]
    }
    return [expr {($dU / ($::BOLTZMANN*[$::thermostatTempCmd])) - $dweight}]
}

# ::xgs::computeNormedWeights
#
# Return the normalized Boltzmann weight of each state at the current
# configuration. These are the conditional probabilities to be used in Gibbs
# sampling procedures.
#
# Arguments:
# ----------
# None
#
# Returns:
# --------
# weights : list of floats
#   A list of normalized weights for each state. 
# 
proc ::xgs::computeNormedWeights {} {
    # Use the max log term as a shift during normalization.
    set logQs [list]
    set logQMax {}
    foreach lambda [getParams] w [getWeights] {
        set U 0.0
        foreach lambdai [getLambdas $lambda]\
                Ui [getAlchEnergies ::energyArray] {
            set U [expr {$U + $lambdai*$Ui}]
        }
        lappend logQs [expr {$w - $U / ($::BOLTZMANN*[$::thermostatTempCmd])}]
        if {[lindex $logQs end] > $logQMax} {
            set logQMax [lindex $logQs end]
        }
    }
    return [normalizeLogWeights $logQs $logQMax]
}

# =============================================================================
# Calibration Routines - necessary to run xgsCalibrate
# =============================================================================
# ::xgs::setupCalibration
#
# Allocate space for the energies that need to be accumulated.
#
proc ::xgs::setupCalibration {} {
    return [lrepeat [getNumStates] [lrepeat 6 0.0]]
}

# ::xgs::accumulateCalibration
#
# Accumulate sums for the current cycle at the current state index.
#
# This procedure should modify the input energyMeans list (use upvar).
#
# Arguments:
# ----------
# energyMeans : list of floats of list of list of floats 
#   List of relevant energy means to be accumulated in each state
# index : int
#   The index of state in which data is being accumulated
#
# Returns:
# --------
# None 
# 
proc ::xgs::accumulateCalibration {energyMeans index} {
    upvar 1 $energyMeans EnergyMeans
    llincr EnergyMeans $index [getAlchEnergies ::energyArray]
    return
}

# ::xgs::computeCalibration
#
# Estimate the state weights - usually by a simple linear response treatment.
#
# Arguments:
# ----------
# energyMeans : list of floats of list of list of floats
#   List of relevant energy means in each state
#
# Returns:
# --------
# weightList : list of floats
#   Estimated weight of each state
#
proc ::xgs::computeCalibration {energyMeans} {
    set weightList [getWeights]
    set numPairs [expr {[getNumStates] - 1}]
    set varList [lrepeat $numPairs 0.0]
    set beta [expr {1.0 / ($::BOLTZMANN*[$::thermostatTempCmd])}]
    for {set i 1} {$i < [getNumStates]} {incr i} {
        set im1 [expr {$i-1}]
        set lambdaim1s [getLambdas [getParams $im1]]
        set lambdais [getLambdas [getParams $i]]
        set Uim1s [lindex $energyMeans $im1]
        set Uis [lindex $energyMeans $i]

        set dUi 0.0
        set dUim1 0.0
        foreach lim1 $lambdaim1s li $lambdais Uim1 $Uim1s Ui $Uis {
             set dl [expr {$li - $lim1}]
             set dUi [expr {$dUi + $dl*$Ui}]
             set dUim1 [expr {$dUim1 + $dl*$Uim1}]
        }
        set wim1 [lindex $weightList $im1]
        lset weightList $i [expr {$wim1 + 0.5*$beta*($dUim1 + $dUi)}]
        lset varList $im1 [expr {$beta*($dUim1 - $dUi)}]
    }
    return [list $weightList $varList]
}

