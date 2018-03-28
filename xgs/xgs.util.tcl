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
# ::xgs::util::optimalWorkStdDev
#
# Return the optimal standard deviation of the instantaneous work for the given
# Gibbs sampling approach. This is optimal in the sense of maximizing the
# round trip rate (or minimizing the round trip time) and thus depends on the
# underlying transition matrix for the given method.
#
# Arguments:
# ----------
# gibbsMethodName : string
#   name of the Gibbs sampling method, any of:
#     stochastic-neighbor - choose a neighbor randomly
#     deterministic-neighbor - choose neighbor in a fixed direction until a
#        a move is rejected
#     convective-neighbor - choose a neighbor in a fixed direction until the
#        endpoint is reached
#     metropolized - choose a new state by independence sampling and then 
#        accept/reject with a modified Metropolis criterion
#
# Returns:
# --------
# sopt : float
#   the standard deviation (in kT units)
#
proc ::xgs::util::optimalWorkStdDev {gibbsMethodName} {
    switch -- [string tolower $gibbsMethodName] {
        "stochastic-neighbor" -
        "metropolized" {
            set sopt 2.381
        }
        "deterministic-neighbor" {
            set sopt 1.729
        }
        "convective-neighbor" {
            set sopt 1.504
        }
        default {
            abort "Unknown Gibbs sampling procedure: $gibbsMethodName"
        }
    }
    return $sopt
}

# ::xgs::util::geometricTemperatures
#
# Return a temperature list with geometric spacing:
#
# T_i = Tmin * (Tmax/Tmin)**[(i - 1) / (M - 1)];    i = 1, ..., M
#
# Arguments:
# ----------
# Tmin : float
#   minimum temperature (in K)
# Tmax : float
#   maximum temperature (in K)
# M : int
#   number of temperatures in ladder (including Tmin and Tmax)
# prec : int (optional, default: 1)
#   round temperatures to "prec" decimal places
#
# Returns:
# --------
# tempList : list
#   a list of the temperatures
#
proc ::xgs::util::geometricTemperatures {Tmin Tmax M {prec 1}} {
    set Ratio [expr {1.0*$Tmax / $Tmin}]
    set Scale [expr {pow(10, int($prec))}]
    set M [expr {int($M)}]
    set tempList [list]
    for {set i 1} {$i <= $M} {incr i} {
        set T [expr {$Tmin*pow($Ratio, [expr ($i - 1.) / ($M - 1.)])}]
        lappend tempList [expr {1.0*round($Scale*$T) / $Scale}]
    }
    return $tempList
}

# ::xgs::util::optimalTempCount
#
# Return the "optimal" (in the sense of minimizing the round trip transit time)
# number of temperatures for a system with the given heat capacity.
#
# This assumes an approximately constant heat capacity, which can only be true 
# over small regions at high temperatures.
#
# Arguments:
# ----------
# heatCapacity : float
#   dimensionless heat capacity (in kB units)
# Tmin : float
#   minimum temperature (in K)
# Tmax : float
#   maximum temperature (in K)
# gibbsMethodName : string
#   name of the Gibbs sampling method (see ::xgs::util::optimalWorkStdDev)
#
# Returns:
# --------
# Mopt : int
#   the number of temperatures
#
proc ::xgs::util::optimalTempCount {heatCapacity Tmin Tmax gibbsMethodName} {
    set sopt [optimalWorkStdDev $gibbsMethodName]
    set s_2 [expr {sqrt($heatCapacity)*log(1.0*$Tmax / $Tmin)}]
    set Mopt [expr {int(round(1 + $s_2/$sopt))}]
    return $Mopt
}

# ::xgs::util::optimalTempLadder 
#
# Return the "optimal" (in the sense of minimizing the round trip transit time)
# temperature ladder for a system with the given (isochoric) heat capacity.
#
# This assumes an approximately constant heat capacity, which can only be true 
# over small regions at high temperatures.
#
# Arguments:
# ----------
# heatCapacity : float
#   dimensionless heat capacity (in kB units)
# Tmin : float
#   minimum temperature (in K)
# Tmax : float
#   maximum temperature (in K)
# gibbsMethodName : string
#   name of the Gibbs sampling method (see ::xgs::util::optimalWorkStdDev)
# prec : int (optional, default: 1)
#   round temperatures to "prec" decimal places
#
# Returns:
# --------
# tempList : list
#   a list of the temperatures
#
proc ::xgs::util::optimalTempLadder {heatCapacity Tmin Tmax gibbsMethodName\
        {prec 1}} {
    set M [optimalTempCount $heatCapacity $Tmin $Tmax $gibbsMethodName]
    return [geometricTemperatures $Tmin $Tmax $M $prec]
}

# ::xgs::util::linearLambdas
#
# Return a lambda list with linear spacing:
#
# l_i = lmin + (lmax - lmin)(i-1)/(M-1)    i = 1, ..., M
#
# Arguments:
# ----------
# lmin : float
#   minimum lambda 
# lmax : float
#   maximum lambda 
# M : int
#   number of lambdas in ladder (including lmin and lmax)
# prec : int (optional, default: 3)
#   round lambda to "prec" decimal places
#
# Returns:
# --------
# lambdaList : list
#   a list of the lambda values
#
proc ::xgs::util::linearLambdas {lmin lmax M {prec 3}} {
    set Scale [expr {pow(10, int($prec))}]
    set M [expr {int($M)}]
    set deltaLambda [expr {($lmax - $lmin) / ($M - 1.)}]
    set lambdaList [list]
    for {set i 1} {$i <= $M} {incr i} {
        set l [expr {$lmin + $deltaLambda*($i - 1)}]
        lappend lambdaList [expr {1.0*round($Scale*$l) / $Scale}]
    }
    return $lambdaList
}

# ::xgs::util::optimalLambdaCount
#
# Return the "optimal" (in the sense of minimizing the round trip transit time)
# number of lambda values for a system with the given alchemical Hessian.
#
# Arguments:
# ----------
# d2fdl2 : float 
#   second derivative of the free energy wrt lambda (in kT units)
# lmin : float
#   minimum lambda 
# lmax : float
#   maximum lambda 
# gibbsMethodName : string
#   name of the Gibbs sampling method (see ::xgs::util::optimalWorkStdDev)
#
# Returns:
# --------
# Mopt : int
#   the number of lambda values
#
proc ::xgs::util::optimalLambdaCount {d2fdl2 lmin lmax gibbsMethodName} {
    set sopt [optimalWorkStdDev $gibbsMethodName]
    set s_2 [expr {sqrt(abs($d2fdl2))}]
    set Mopt [expr {int(round(1 + $s_2/$sopt))}]
    return $Mopt 
}

# ::xgs::util::optimalLambdaLadder 
#
# Return the "optimal" (in the sense of minimizing the round trip transit time)
# lambda ladder for a system with the given alchemical Hessian.
#
# Arguments:
# ----------
# d2fdl2 : float 
#   second derivative of the free energy wrt lambda (in kT units)
# lmin : float
#   minimum lambda 
# lmax : float
#   maximum lambda 
# gibbsMethodName : string
#   name of the Gibbs sampling method (see ::xgs::util::optimalWorkStdDev)
# prec : int (optional, default: 3)
#   round lambda to "prec" decimal places
#
# Returns:
# --------
# lambdaList : list
#   a list of the lambda values
#
proc ::xgs::util::optimalLambdaLadder {d2fdl2 lmin lmax gibbsMethodName\
        {prec 3}} {
    set M [optimalLambdaCount $d2fdl2 $lmin $lmax $gibbsMethodName]
    return [linearLambdas $lmin $lmax $M $prec]
}

