# xgs.core.tcl
#
# This file provides the core functionality for eXpanded ensemble simulations 
# with Gibbs Sampling (XGS) in NAMD. It is not meant to be used directly for
# simulations, but rather to be expanded upon in order to implement the 
# specific type of parameter sampling that is desired.
#
package require Tcl 8.5

source [file join [file dirname [info script]] "namdtcl.tcl"]
source [file join [file dirname [info script]] "numtcl.tcl"]
source [file join [file dirname [info script]] "namdmcmc.tcl"]

namespace eval ::xgs {
    variable TITLE "XGS)"
    variable TYPE ""
    #   States are defined by a parameter (or parameters) stored in ParamLadder
    # and a weight stored in StateWeights. Each state then has a unique 
    # identifier equal to the index in these lists. For dimensions greater than
    # one an orthorhombic space is assumed and the number of parameters along
    # each axis is stored in ParamDims. The number of states is thus equal to
    # the product of all elements in ParamDims.
    #
    variable ParamLadder
    variable ParamDims
    variable StateWeights [list]
    # Because the initial index can be inferred in multiple ways, the initial
    # and current state index are stored separately (see xgsSetup).
    variable initialStateIndex -1
    variable stateIndex -1
    # For history dependent neighbor schemes
    variable stepIncr
    variable prevAxis -1
    # Type of Gibbs sampling protocol to use.
    variable GibbsMethodName "metropolized"
    variable GibbsSamplingCmd
    # Information for reading/writing restart files.
    variable restartFilename ""
    variable restartFreq 0

    namespace export xgsRun xgsCalibrate\
            xgsParamLadder xgsStateWeights xgsGibbsMethod\
            xgsRestartFile xgsRestartFreq
}

# Check that a proc is defined - this is used to make sure that a specific
# method has implemented the procs that XGS requires.
#
proc ::xgs::procExists {procName} {
    return [expr {[llength [info procs $procName]] > 0}]
}

# =============================================================================
# Main run commands
# =============================================================================
# ::xgs::xgsRun
#
# Run a given number of MD/MC cycles using a fixed number of steps per cycle.
#
# In order to get the most sensible output from NAMD, the MC step is executed
# first in each cycle, followed by MD.
# 
# Arguments:
# ----------
# numsteps : int
#   The number of MD steps per cycle
# numcycles : int
#   The number of MC cycles
#
# Returns:
# --------
# None
#
proc ::xgs::xgsRun {numsteps numcycles} {
    xgsSetup
    run 0
    set acceptCount 0
    set firstCycle 1
    set lastCycle [expr {$firstCycle + int($numcycles) - 1}]
    for {set cycle $firstCycle} {$cycle <= $lastCycle} {incr cycle} {
        # (1) Perform a MC move in parameter space.
        lassign [sampleNewState] index accept
        incr acceptCount $accept
        set param [getParams $index]
        xgsPrint "cycle $cycle : parameter(s) $param : index $index"
        # (2) Perform MD with the current parameters.
        run norepeat [expr int($numsteps)]
        writeRestart "[outputname].xgsrst" $cycle
    }
    writeRestart force "[outputname].xgsrst" $cycle
    set acceptanceRate [expr {100.*$acceptCount / $numcycles}]
    xgsPrint [format "acceptance rate = %5.1f %%" $acceptanceRate]
    return
}

# ::xgs::xgsCalibrate
#
# Scan through all of the states in order to make a quick estimate of the state
# weights. Data is collected at intervals of "outputEnergies", so the number of
# steps per state must be a whole number multiple of this.
#
# Upon completion, an XGS restart file is written that contains all of the 
# relevant calibration data.
# 
# Arguments:
# ----------
# numSteps : int
#   The number of MD steps to collect data for in each state
# numEquilSteps : int (optional, default: 0)
#   The number of MD steps to equilibrate in each state
#
# Returns:
# --------
# weights : list
#   The estimated free energy weights for each state
# variances : list
#   The estimated work variances _between_ neighboring pairs of states
#
proc ::xgs::xgsCalibrate {numsteps {numEquilSteps 0}} {
    # Check that the necessary subroutines are defined.
    foreach procName [list setupCalibration accumulateCalibration\
                           computeCalibration] {
        if {![procExists $procName]} {
            xgsAbort "$procName is not defined!"
        }
    }
    xgsStateIndex 0
    xgsSetup
    if {[expr {$numsteps % [outputEnergies]}] != 0} {
        xgsAbort "Number of calibration steps must be a whole number multiple"\
                 "of outputEnergies"
    }
    set numsamples [expr {int($numsteps / [outputEnergies])}]
    # (1) Determine mean energies in each state.
    set energyMeans [setupCalibration]
    set index 0
    foreach paramList [getParams] {
        if {$numEquilSteps} {
            xgsPrint "Equilibrating at $paramList"
            if {$index} {
                updateState [expr {$index - 1}] $index
                firstTimestep 0
            }
            run $numEquilSteps
        }
        xgsPrint "Calibrating at $paramList"
        firstTimestep 0
        run 0
        for {set n 0} {$n < $numsamples} {incr n} {
            run norepeat [outputEnergies]
            storeEnergies
            accumulateCalibration energyMeans $index
        }
        lmultiply energyMeans $index [expr {1.0 / $numsamples}]
        incr index
    }
    # (2) Compute the state weights.
    lassign [computeCalibration $energyMeans] weights variances
    xgsStateWeights $weights
    # (3) Report to the user - mimic the energy logging format.
    xgsPrint "Calibration Report:"
    xgsPrint "Parameters:"
    xgsPrintEnergyLabels [getParams]
    xgsPrint "Weights:"
    xgsPrintEnergies $weights
    xgsPrint "Variances:"
    xgsPrintEnergies $variances
    xgsPrint "Results written to XGS restart [outputname].xgsrst"
    writeRestart force "[outputname].xgsrst" 0
    return [list $weights $variances $energyMeans]
}

# =============================================================================
# Setter functions to be used as NAMD keywords
# =============================================================================
# ::xgs::xgsGibbsMethod
#
# Set the desired Gibbs sampling method. Currently supported options are:
#
# neighbor - nearest neighbor sampling (stochastic selection of parameter axis 
#              and step direction with reflective boundaries)
#
# metropolized* - Metropolized independence sampling
#
# * - default
# 
proc ::xgs::xgsGibbsMethod {methodName} {
    set ::xgs::GibbsMethodName [string trim $methodName]
    return
}

# ::xgs::xgsParamLadder
#
# Set or reset the state parameter ladder. 
#
# These are NOT sorted or processed in any way and will be used _as is_ for 
# determining neighbors.
# 
proc ::xgs::xgsParamLadder {args} {
    variable ::xgs::ParamLadder [list]
    variable ::xgs::ParamDims
    variable ::xgs::stepIncr
    set dim [llength $args]
    if {$dim == 1} {
        set ParamDims 0
        foreach param {*}$args {
            checkIsNumeric xgsParamLadder $param
            lappend ParamLadder $param
            incr ParamDims
        }
    } elseif {$dim == 2} {
        set ParamDims [list]
        foreach paramList $args {
            lappend ParamDims [llength $paramList]
        }
        foreach param1 [lindex $args 1] {
            checkIsNumeric xgsParamLadder $param1
            foreach param0 [lindex $args 0] {
                checkIsNumeric xgsParamLadder $param0
                lappend ParamLadder [list $param0 $param1]
            }
        }
    } else {
        xgsAbort "Parameter spaces higher than 2 dimensions are not currently"\
                 "supported."
    }
    set stepIncr [lrepeat $dim 1]
    return
}

# ::xgs::xgsStateWeights
#
# Set or reset the state weights (in kT units).
#
# Using the exact free energies should lead to a uniform distribution across 
# all states.
# 
proc ::xgs::xgsStateWeights {weights} {
    variable ::xgs::StateWeights [list]
    foreach weight $weights {
        checkIsNumeric xgsStateWeights $weight
        lappend StateWeights $weight
    }
    return
}

# ::xgs::xgsStateIndex
#
# Set the initial state index.
#
proc ::xgs::xgsStateIndex {initialIndex} {
    checkIsNumeric xgsStateIndex $initialIndex
    # Note that we do not use the standard setter here, because that requires
    # that the state space be defined - save that check for later.
    variable ::xgs::initialStateIndex [expr {int($initialIndex)}]
    return
}

# ::xgs::xgsRestartFile
#
# Initialize or re-initialize settings from an XGS restart file.
#
proc ::xgs::xgsRestartFile {filename} {
    variable ::xgs::restartFilename $filename
    return
}

# ::xgs::xgsRestartFreq
#
# Set the frequency (in MC cycles) at which restart files are written. The
# default, 0, only writes a restart upon normal completion.
#
proc ::xgs::xgsRestartFreq {frequency} {
    checkIsNotNegative xgsRestartFreq frequency
    variable ::xgs::restartFreq [expr {int($frequency)}]
    return
}

# =============================================================================
# Convenience procedures for I/O
# =============================================================================
# print statement with a pretty title
proc ::xgs::xgsPrint {args} {
    print "$::xgs::TITLE [join $args " "]"
}

# abort statement with a pretty title
proc ::xgs::xgsAbort {args} {
    abort "$::xgs::TITLE [join $args " "]"
}

# Print a list of energy labels - use the standard NAMD print format.
#
proc ::xgs::xgsPrintEnergyLabels {energyLabelList} {
    set str "     "
    set count 1
    foreach energyLabel $energyLabelList {
        set str [format "%s %14s" $str $energyLabel]
        incr count
        if {$count == 5} {
            set count 0
            set str "$str     "
        }
    }
    xgsPrint $str
}

# Print a list of energy values - use the standard NAMD print format.
#
proc ::xgs::xgsPrintEnergies {energyList} {
    set str "     "
    set count 1
    foreach energy $energyList {
        set str [format "%s % 14.4f" $str $energy]
        incr count
        if {$count == 5} {
            set count 0
            set str "$str     "
        }
    }
    xgsPrint $str
}

# ::xgs::writeRestart
#
# Write the current state information as well as the parameter ladder and state
# weights to file based on the restart frequency. If the "force" keyword 
# precedes the arguments, ignore the restart frequency and write no matter 
# what.
#
# Arguments:
# ----------
# restartFilename : string
#   Name of (JSON format) restart file to write
# cycle : int
#   Last neMD/MC cycle index
#
# Returns:
# --------
# true if a restart file was written, else false
#
proc ::xgs::writeRestart {args} {
    variable ::xgs::restartFreq
    if {[string match [lindex $args 0] force]} {
        set restartFilename [lindex $args 1]
        set cycle [expr {int([lindex $args 2])}]
    } else {
        set restartFilename [lindex $args 0]
        set cycle [expr {int([lindex $args 1])}]
        if {!($restartFreq) || ![expr {$cycle % $restartFreq}]} {
            return 0
        }
    }
    namdFileBackup $restartFilename
    set RestartFile [open $restartFilename "w"]
    puts $RestartFile "xgsStateIndex [getStateIndex]"
    for {set axis 0} {$axis < [getNumDims]} {incr axis} {
        puts $RestartFile "xgsParamLadder [getParamAxis $axis]"
    }
    puts $RestartFile "xgsStateWeights [getWeights]"
    close $RestartFile
    return 1
}

# ::xgs::readRestart
#
# Read in stored state information from a restart file. This includes all
# settings that must remain fixed during a simulation in order to yield a
# proper ensemble. That is, if any of these parameters are changed manually,
# then this is no longer a valid continuation of a previous simulation.
#
proc ::xgs::readRestart {restartFilename} {
    set RestartFile [open $restartFilename "r"]
    set Restart [string trim [read $RestartFile]]
    close $RestartFile
    set ParamLadderArgs [list]
    foreach line [split $Restart "\n"] {
        set tokens [split $line]
        set key [lindex $tokens 0]
        set val [join [lrange $tokens 1 end]]
        switch -- $key {
            "xgsStateIndex" {
                xgsStateIndex $val 
            }
            "xgsParamLadder" {
                # This is a bit of a hack to account for the fact that lists
                # and lists of lists are not easily differentiated...
                lappend ParamLadderArgs $val
            }
            "xgsStateWeights" {
                xgsStateWeights $val
            }
            default {
                xgsPrint "WARNING: Unknown restart keyword: $key"
            }
        }
    }
    if {[llength $ParamLadderArgs] > 0} {
        xgsParamLadder {*}$ParamLadderArgs
    }    
    return
}

# =============================================================================
# Gibbs Sampling Routines 
# =============================================================================
# ::xgs::stochasticNeighborSampling
#
# Propose and accept/reject a change to a neighboring state. The candidate
# state is chosen randomly amongst neighbors.
#
# Arguments:
# ----------
# oldIndex : int
#   The index of the current state in the lambda ladder
#
# Returns:
# --------
# newIndex : int
#   The new state index after the Metropolis accept/reject step
# 
proc ::xgs::stochasticNeighborSampling {oldIndex} {
    set step [expr {int(pow(-1, randint(0, 1)))}]
    set axis [expr {randint(0, [expr {[getNumDims] - 1}])}]
    return [neighborSampling $oldIndex $step $axis]
}

# ::xgs::deterministicNeighborSampling
#
# Propose and accept/reject a change to a neighboring state. The candidate is
# chosen by a deterministic up/down sequence where the direction changes when
# a proposal is rejected.
#
# Arguments:
# ----------
# oldIndex : int
#   The index of the current state in the lambda ladder
#
# Returns:
# --------
# newIndex : int
#   The new state index after the Metropolis accept/reject step
# 
proc ::xgs::deterministicNeighborSampling {oldIndex} {
    variable ::xgs::prevAxis
    incr prevAxis
    set prevAxis [expr {($prevAxis == [getNumDims]) ? 0 : $prevAxis}]
    set axis $prevAxis
    set step [lindex $::xgs::stepIncr $axis]
    return [neighborSampling $oldIndex $step $axis]
}

# ::xgs::convectiveNeighborSampling
#
# Propose and accept/reject a change to a neighboring state. The candidate is
# always driven towards the next endpoint.
#
# Arguments:
# ----------
# oldIndex : int
#   The index of the current state in the lambda ladder
#
# Returns:
# --------
# newIndex : int
#   The new state index after the Metropolis accept/reject step
# 
proc ::xgs::convectiveNeighborSampling {oldIndex} {
    return [deterministicNeighborSampling $oldIndex]
}

# ::xgs::neighborSampling
#
# Propose and accept/reject a change to a neighboring state.
#
# Arguments:
# ----------
# oldIndex : int
#   The index of the current state in the lambda ladder
# step : int
#   The step size and direction of the new neighbor
# axis : int
#   The (if applicable) parameter axis to find a neighbor on
#
# Returns:
# --------
# newIndex : int
#   The new state index after the Metropolis accept/reject step
# 
proc ::xgs::neighborSampling {oldIndex step axis} {
    # For now, 1 and 2D are treated separately - it might be fairly easy to
    # generalize to d dimensions, but not today...
    set du 0.0
    switch -- [getNumDims] {
        1 {
            set newIndex [expr {$oldIndex + $step}]
            set tmpIndex [correctForBoundary $newIndex [getNumStates]]
            if {$newIndex != $tmpIndex} {
                set newIndex $tmpIndex
                set du $::LN2
            }
        }
        2 {
            lassign [getParamDim] N0 N1
            set oldi [expr {$oldIndex % $N0}]
            set oldj [expr {$oldi / $N0}]
            if {$axis == 0} {
                set newi [expr {$oldi + $step}]
                set tmpi [correctForBoundary $newi $N0]
                if {$newi != $tmpi} {
                    set newi $tmpi
                    set du $::LN2
                }
                set newIndex [expr {$oldj*$N0 + $newi}]
            } elseif {$axis == 1} {
                set newj [expr {$oldj + $step}]
                set tmpj [correctForBoundary $newj $N1]
                if {$newj != $tmpj} {
                    set newj $tmpj
                    set du $::LN2
                }
                set newIndex [expr {$newj*$N0 + $oldi}]
            } else {
                xgsAbort "Error when choosing a parameter axis"
            }
        }
        default {
            xgsAbort "Bad parameter dimensions! [getNumDims], expected 1 or 2"
        }
    }
    storeEnergies
    set du [expr {$du + [computePairPotential $oldIndex $newIndex]}]
    return [reportAndUpdate $du $oldIndex $newIndex]
}

# Convenient shorthand for reflective boundary conditions.
#
# Briefly, always move away from a boundary. However, since this is a 100/0
# proposal instead of 50/50, the Metropolis criterion has to be corrected by a
# factor or 1/2. This is equivalent to an energy shift of ln(2) kT units.
#
proc ::xgs::correctForBoundary {index dim} {
    if {$index < 0} {
        return 1
    } elseif {$index >= $dim} {
        return [expr {$dim - 2}]
    } else {
        return $index
    }
}

# ::xgs::metropolizedSampling
#
# Propose and accept/reject a new index using the Metropolized independence
# criteria. Return the index of the new state. Note that the current state is
# _never_ proposed, but the proposed state can be rejected (unlike 
# conventional independence sampling).
#
# Arguments:
# ----------
# oldIndex : int
#   The index of the current state in the lambda ladder
# tol : float (default: 1e-3)
#   Tolerance for ignoring the proposal step. If the _combined_ probability of
#   all other states is less than tol then the current state is just accepted.
#
# Returns:
# --------
# newIndex : int
#   The new state index after the Metropolis accept/reject step
#
proc ::xgs::metropolizedSampling {oldIndex {tol 1e-3}} {
    # Save the probability of the current state, i as p_i. The proposal
    # probability for state i is then set to 0 and, for any other state j, it 
    # is set to:
    #
    # p_j / (1 - p_i)
    #
    # Since the weights are already normalized, we can save computations by
    # zeroing out p_i and noting that the normalization of everything else is 
    # just (1 - p_i).
    #
    # For stability (and to save some time), don't play this game if the
    # aggregate probability of all other states is less than the tolerance.
    #
    storeEnergies
    set Weights [computeNormedWeights]
    set pic [expr {1. - [lindex $Weights $oldIndex]}] ;# the "complement" of pi
    if {$pic < $tol} {
        set newIndex $oldIndex
        set du 0.0
    } else {
        lset Weights $oldIndex 0.0
        lassign [choice $Weights $Weights] pj newIndex
        set du [expr {log((1. - $pj) / $pic)}]
    }
    return [reportAndUpdate $du $oldIndex $newIndex]
}

# ::xgs::sampleNewState
#
# Select a new state using the selected Gibbs sampling method.
#
proc ::xgs::sampleNewState {} {
    return [$::xgs::GibbsSamplingCmd [getStateIndex]]
}

# ::xgs::xgsSetup
#
# Setup parameters for Gibbs sampling. Also do basic sanity checks.
#
# Arguments:
# ----------
# None
#
# Returns:
# --------
# None
#
proc ::xgs::xgsSetup {} {
    variable ::xgs::TYPE
    variable ::xgs::ParamLadder
    variable ::xgs::StateWeights
    variable ::xgs::GibbsMethodName
    variable ::xgs::GibbsSamplingCmd
    variable ::xgs::restartFilename
    variable ::xgs::initialStateIndex

    if {[string length $restartFilename] > 0} {
        xgsPrint "Reading XGS settings from $restartFilename"
        readRestart $restartFilename
    }
    # Sanity check the parameter ladder and state weights:
    #   1) Default to all zero weights if none were given.
    #   2) Check that the same number of states are indicated.
    #
    if {[llength $StateWeights] == 0} {
        foreach param $ParamLadder {
            lappend StateWeights 0.0
        }
    }
    if {[string length $TYPE] > 0} {
        xgsPrint "simulation type: $TYPE"
    }
    xgsPrint "state parameters: $ParamLadder"
    xgsPrint "state weights: $StateWeights"
    if {[llength $StateWeights] != [llength $ParamLadder]} {
        xgsAbort "Parameter ladder and weight specifications do not match!"
    }
    # NB setStateIndex requires ParamLadder to be defined, so we can't do this
    # until this point.
    #
    if {$initialStateIndex != -1} {
        setStateIndex $initialStateIndex
    } else {
        # Default to state 0.
        setStateIndex 0
    }
    # Select the Gibbs sampling routine.
    #
    switch -- [string tolower $GibbsMethodName] {
        "stochastic-neighbor" {
            xgsPrint "Using stochastic nearest neighbor sampling"
            set GibbsSamplingCmd stochasticNeighborSampling
            set GibbsWeightCmd computePairPotential
        }
        "deterministic-neighbor" {
            xgsPrint "Using deterministic nearest neighbor sampling"
            set GibbsSamplingCmd deterministicNeighborSampling
            set GibbsWeightCmd computePairPotential
        }
        "convective-neighbor" {
            xgsPrint "Using convective nearest neighbor sampling"
            set GibbsSamplingCmd convectiveNeighborSampling
            set GibbsWeightCmd computePairPotential
        }
        "metropolized" {
            xgsPrint "Using Metropolized independence sampling"
            set GibbsSamplingCmd metropolizedSampling
            set GibbsWeightCmd computeNormedWeights
        }
        default {
            xgsAbort "Unknown Gibbs sampling procedure: $GibbsMethodName"
        }
    }
    if {![procExists $GibbsWeightCmd]} {
        xgsAbort "Gibbs method $GibbsSamplingCmd selected, but"\
                 "$GibbsWeightCmd is not defined"
    }
    if {![procExists parameterSetup] || ![procExists updateState]} {
        xgsAbort "parameterSetup and updateState are not defined!"
    }
    # Select the thermostat temperature to be used for MC.
    #
    getThermostat
    if {!$::thermostatIsSet} { 
       xgsAbort "A thermostat is required"
    }
    xgsPrint "Using $::thermostatName thermostat"
    getBarostat
    if {$::barostatIsSet} {
        xgsPrint "Using $::barostatName barostat"
    }

    # Finish with method specific setup.
    #
    parameterSetup
    callback energyCallback
    return
}

# ::xgs::reportAndUpdate
#
# Accept/reject newIndex in place of oldIndex using the Metropolis criterion on
# the reduced energy du. Return the index that was chosen and appropriately
# update the simulation parameters.
#
# Arguments:
# ----------
# du : float
#   The reduced energy used in the Monte Carlo criteria
# oldIndex : int
#   The index of the old state in the parameter ladder
# newIndex : int
#   The index of the new state in the parameter ladder
#
# Returns:
# --------
# currentIndex : int
#   The current index (will be either oldIndex or newIndex)
#
proc ::xgs::reportAndUpdate {du oldIndex newIndex} {
    variable ::xgs::GibbsMethodName
    variable ::xgs::stepIncr
    variable ::xgs::prevAxis

    if {$oldIndex != $newIndex} {
        set oldParam [getParams $oldIndex]
        set newParam [getParams $newIndex]
        set accept [metropolisAcceptance $du]
        xgsPrint [format "exchange %s --> %s : dE = %.4f kT :\
                          accept = %d" $oldParam $newParam $du $accept
              ]
    } else {
        set accept 0
        xgsPrint "no exchange"
    }

    set prevIncr [lindex $stepIncr $prevAxis]
    if {$accept} {
        updateState $oldIndex $newIndex
        setStateIndex $newIndex
        # This is only needed for history dependent neighbor sampling
        if {[string match -nocase $GibbsMethodName "convective-neighbor"]} {
            # prevIncr is 1 if we are trying to hit the upper endpoint and -1
            # if we are trying to hit the lower endpoint on this axis.
            lassign [getParamDim] N0
            set hit 0
            if {$prevAxis == 0} {
                if {$prevIncr == 1 && ![expr {($newIndex+1)%$N0}]} {
                    set hit 1
                } elseif {$prevIncr == -1 && ![expr {$newIndex%$N0}]} {
                    set hit 1
                }
            } elseif {$prevAxis == 1} {
                # WARNING! THIS IS ONLY CORRECT FOR 2 OR FEWER DIMENSIONS!
                if {$prevIncr == 1 
                    && [expr {$newIndex + $N0}] >= [getNumStates]} {
                    set hit 1
                } elseif {$prevIncr == -1 && $newIndex < $N0} { 
                    set hit 1
                }
            }
            if {$hit} {
                lset stepIncr $prevAxis [expr {-1*$prevIncr}]
            }
        }
    } else {
        # This is only needed for history dependent neighbor sampling
        if {[string match -nocase $GibbsMethodName "deterministic-neighbor"]} { 
            lset stepIncr $prevAxis [expr {-1*$prevIncr}]
        }
    }
    return [list [getStateIndex] $accept]
}

# =============================================================================
# Setter functions
# =============================================================================
# ::xgs::setStateIndex
#
# Set the current state index - THIS DOES NOT UPDATE THE STATE PARAMETERS!
#
proc ::xgs::setStateIndex {index} {
    if {$index < 0 || [getNumStates] <= $index} {
        xgsAbort "Invalid state index $index, must be in range"\
                 "\[0, [getNumStates]\)."
    }
    variable ::xgs::stateIndex $index
    return
}

# ::xgs::setTemp
#
# Set/reset the thermostat temperature(s).
#
proc ::xgs::setTemp {newTemp} {
    $::thermostatTempCmd $newTemp
    if {$::barostatIsSet && [$::barostatCmd]} {
        if {[string length $::barostatTempCmd] > 0} {
            $::barostatTempCmd $newTemp
        }
    }
    return
}

# =============================================================================
# Getter functions
#
# These are largely unnecessary, but cut down on "variable" statements.
# =============================================================================
# ::xgs::getStateIndex
#
# Return the current state index.
#
proc ::xgs::getStateIndex {} {
    return $::xgs::stateIndex
}

# ::xgs::getParams
#
# Return the parameter(s) of a state given its index. If no index is given,
# return a list of all of the parameters.
#
proc ::xgs::getParams {args} {
    variable ::xgs::ParamLadder
    if {[llength $args] == 0} {
        return $ParamLadder
    } else {
        return [lindex $ParamLadder {*}$args]
    }
}

# ::xgs::getParamAxis
#
# Return the parameters along a given axis.
#
proc ::xgs::getParamAxis {axis} {
    set ParamAxis [list]
    switch -- [getNumDims] {
        1 {
            if {$axis == 0} {
                set ParamAxis [getParams]
            } else {
                xgsAbort "Bad parameter axis $axis."
            }
        }
        2 {
            if {$axis == 0} {
                for {set i 0} {$i < [getParamDim 0]} {incr i} {
                    lappend ParamAxis [getParams $i 0]
                }
            } elseif {$axis == 1} {
                for {set j 0} {$j < [getParamDim 1]} {incr j} {
                    set ij [expr {$j*[getParamDim 0]}]
                    lappend ParamAxis [getParams $ij 1]
                }
            } else {
                xgsAbort "Bad parameter axis $axis."
            }
        }
        default {
            xgsAbort "Parameter spaces higher than 2 dimensions are not"\
                     "currently supported."
        }
    }
    return $ParamAxis
}

# ::xgs::getParamIndex
#
# Return the index of a state given its parameter(s).
#
proc ::xgs::getParamIndex {params} {
    return [lsearch $::xgs::ParamLadder $params]
}

# ::xgs::getWeights
#
# Return the weight of a state given its index. If no index is given, return a
# list of all of the weights.
#
proc ::xgs::getWeights {args} {
    variable ::xgs::StateWeights
    if {[llength $args] == 0} {
        return $StateWeights
    } else {
        return [lindex $StateWeights {*}$args]
    }
}

# ::xgs::getNumStates
#
# Return the number of a states.
#
proc ::xgs::getNumStates {} {
    return [llength $::xgs::ParamLadder]
}

# ::xgs::getParamDim
#
# Return the number of parameters along an axis.
#
proc ::xgs::getParamDim {args} {
    variable ::xgs::ParamDims
    if {[llength $args] == 0} {
        return $::xgs::ParamDims
    } else {
        return [lindex $::xgs::ParamDims {*}$args]
    }
}

# ::xgs::getNumDims
#
# Return the number of parameter axes.
#
proc ::xgs::getNumDims {} {
    return [llength $::xgs::ParamDims]
}

