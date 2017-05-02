# xgs.st.tcl
#
# This file implements XGS protocols for running simulated tempering (ST).
# Rather than require NAMD users to handle the ::xgs namespace directly, only
# those parts that are relevant to the user are imported/exported in st.tcl. 
#
package require Tcl 8.5

source [file join [file dirname [info script]] "xgs.core.tcl"]
source [file join [file dirname [info script]] "xgs.util.tcl"]

namespace eval ::xgs {
    variable TYPE ST
    namespace import ::xgs::util::geometricTemperatures \
            ::xgs::util::optimalTempCount ::xgs::util::optimalTempLadder
    namespace export stCalibrate geometricTemperatures \
            optimalTempCount optimalTempLadder
}

# Here we equate the instantaneous work fluctutations needed to change the
# temperature with the heat capacity (which is a constant between temperatures
# within linear response). Further assuming that this is constant over _all_
# temperatures leads to a single estimate. This should converge much faster
# than data from two temperatures, but still probably not too well. In any
# event, even an estimate with 10% error should be pretty good for determining
# an optimal temperature ladder (an error of roughly +/- 2 in the ladder size).
#
proc ::xgs::stCalibrate {numsteps {numEquilSteps 0}} {
    lassign [xgsCalibrate $numsteps $numEquilSteps] weightList varList
    set tempList [getParams]
    set TiList [lrange $tempList 0 end-1]
    set Tip1List [lrange $tempList 1 end]
    set c 0.0
    foreach var $varList Ti $TiList Tip1 $Tip1List {
        set alp [expr {$Tip1/$Ti}]
        set c [expr {$c + $alp*$var/($alp-1)**2}]
    }
    set c [expr {$c / [llength $varList]}]
    set Tmin [lindex $tempList 0]
    set Tmax [lindex $tempList end]
    set N [optimalTempCount $c $Tmin $Tmax]
    set tladder [optimalTempLadder $c $Tmin $Tmax]
    xgsPrint "Properties from linear response:"
    xgsPrint [format "excess heat capacity (kB units): %.2e" $c]
    xgsPrint "optimal temperature count: $N"
    xgsPrint "optimal temperature ladder: $tladder"
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
    set Temp [getParams [getStateIndex]]
    setTemp $Temp
    setStateIndex [getParamIndex $Temp] 
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
    set oldTemp [getParams $oldIndex]
    set newTemp [getParams $newIndex]
    rescalevels [expr {sqrt($newTemp/$oldTemp)}]
    setTemp $newTemp
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
    set oldTemp [getParams $oldIndex]
    set newTemp [getParams $newIndex]
    set dT [expr {$oldTemp - $newTemp}]
    set dbeta [expr {$dT / ($::BOLTZMANN*$oldTemp*$newTemp)}]
    set U $::energyArray(POTENTIAL)
    if {[constPresOn]} {
        set P [$::barostatPresCmd]
        set V $::energyArray(VOLUME)
        set U [expr {$U + $P*$V / $::PRESSUREFACTOR}]
    }
    set dweight [expr {[getWeights $newIndex] - [getWeights $oldIndex]}]
    return [expr {$dbeta*$U - $dweight}]
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
    set U $::energyArray(POTENTIAL)
    if {[constPresOn]} {
        set P [$::barostatPresCmd]  
        set V $::energyArray(VOLUME)
        set U [expr {$U + $P*$V / $::PRESSUREFACTOR}]
    }
    set logQs [list]
    set logQMax {}
    foreach T [getParams] w [getWeights] {
        lappend logQs [expr {$w - $U/($::BOLTZMANN*$T)}]
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
    return [lrepeat [getNumStates] 0.0]
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
    set U $::energyArray(POTENTIAL)
    if {[constPresOn]} {
        set P [$::barostatPresCmd]
        set V $::energyArray(VOLUME)
        set U [expr {$U + $P*$V / $::PRESSUREFACTOR}]
    } 
    lincr EnergyMeans $index $U
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
    for {set i 1} {$i < [getNumStates]} {incr i} {
        set im1 [expr {$i-1}]
        set Tim1 [getParams $im1]
        set Eim1 [lindex $energyMeans $im1 0]
        set Ti [getParams $i]
        set Ei [lindex $energyMeans $i 0]
        set dbeta [expr {($Tim1 - $Ti) / ($::BOLTZMANN*$Tim1*$Ti)}]
        set wim1 [lindex $weightList $im1]
        lset weightList $i [expr {$wim1 + 0.5*$dbeta*($Ei + $Eim1)}]
        lset varList $im1 [expr {$dbeta*($Eim1 - $Ei)}]
    }
    return [list $weightList $varList]
}

