# xgs.util.tcl
#
# This file provides additional utilities for expanded ensemble simulations in
# NAMD. These might be applicable for multiple simulation types, but not for
# all simulation types. As such these are less general than the procedures in
# xgs.core.tcl but more general than what might be found in xgs.st.tcl, for
# example.
#
namespace eval ::xgs::util {
    namespace export * 
}

# =============================================================================
# Utilities for allocating parameters and states
# =============================================================================
# geometricTemperatures
#
# Return a temperature list with geometric spacing:
#
# T_i = Tmin * (Tmax/Tmin)**[(i - 1) / (N - 1)];    i = 1, ..., N
#
# Arguments:
# ----------
# Tmin : float
#   minimum temperature (in K)
# Tmax : float
#   maximum temperature (in K)
# N : int
#   number of temperatures in ladder (including Tmin and Tmax)
# prec : int (optional, default: 1)
#   round temperatures to "prec" decimal places
#
# Returns:
# --------
# tempList : list
#   a list of the temperatures
#
proc ::xgs::util::geometricTemperatures {Tmin Tmax N {prec 1}} {
    set tempList [list]
    set Ratio [expr {1.0*$Tmax / $Tmin}]
    set Scale [expr {pow(10, int($prec))}]
    set N [expr {int($N)}]
    for {set i 1} {$i <= $N} {incr i} {
        set T [expr {$Tmin*pow($Ratio, [expr ($i - 1.) / ($N - 1.)])}]
        lappend tempList [expr {1.0*round($Scale*$T) / $Scale}]
    }
    return $tempList
}

# optimalTempLadder 
#
# Return the "optimal" (in the sense of minimizing the round trip transit time)
# temperature ladder for a system with the given (isochoric) heat capacity.
#
# This assumes an approximately constant heat capacity, which can only be true 
# over small regions at high temperatures.
#
# Reference: 
# (1) R. Denschlag, M. Lingenheil and P. Tavan "Optimal Temperature Ladders in
#   Replica Exchange Simulations" Chem. Phys. Lett. 2009, 472, 193.
#
# Arguments:
# ----------
# heatCapacity : float
#   dimensionless heat capacity (in kB units)
# Tmin : float
#   minimum temperature (in K)
# Tmax : float
#   maximum temperature (in K)
# prec : int (optional, default: 1)
#   round temperatures to "prec" decimal places
# linPenalty : bool (optional, default: true)
#   if true, apply a linear penalty in the ladder size - this is meant to
#   balance between a high transition rate and the cost of sampling more states
#
# Returns:
# --------
# tempList : list
#   a list of the temperatures
#
proc ::xgs::util::optimalTempLadder {heatCapacity Tmin Tmax {prec 1}
        {linPenalty true}} {
    if {$linPenalty} {
        set factor [expr {0.594*sqrt(0.5*$heatCapacity) - 0.5}]
    } else {
        set factor [expr {sqrt(0.5*$heatCapacity) / 1.068 - 0.5}]
    }
    set N [expr {int(round(1 + $factor*log($Tmax / $Tmin)))}]
    return [geometricTemperatures $Tmin $Tmax $N $prec]
}

