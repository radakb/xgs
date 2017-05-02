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
# ::xgs::util::geometricTemperatures
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

# ::xgs::util::optimalTempCount
#
# Return the "optimal" (in the sense of minimizing the round trip transit time)
# number of temperatures for a system with the given (isochoric) heat capacity.
#
# This assumes an approximately constant heat capacity, which can only be true 
# over small regions at high temperatures.
#
# Reference:
# (1) W. Nadler and U. Hansmann "Optimized Explicit-Solvent Replica Exchange
#   Molecular Dynamics from Scratch" J. Phys. Chem. B 2008, 112, 10386.
# (2) R. Denschlag, M. Lingenheil and P. Tavan "Optimal Temperature Ladders in
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
#
# Returns:
# --------
# Nopt : int
#   the number of temperatures
#
proc ::xgs::util::optimalTempCount {heatCapacity Tmin Tmax} {
    set factor [expr {0.594*sqrt(0.5*$heatCapacity)}]
    set N [expr {int(round(1 + $factor*log($Tmax / $Tmin)))}]
    return $N
}

# ::xgs::util::optimalTempLadder 
#
# Return the "optimal" (in the sense of minimizing the round trip transit time)
# temperature ladder for a system with the given (isochoric) heat capacity.
#
# This assumes an approximately constant heat capacity, which can only be true 
# over small regions at high temperatures.
#
# Reference:
# (1) W. Nadler and U. Hansmann "Optimized Explicit-Solvent Replica Exchange
#   Molecular Dynamics from Scratch" J. Phys. Chem. B 2008, 112, 10386.
# (2) R. Denschlag, M. Lingenheil and P. Tavan "Optimal Temperature Ladders in
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
#
# Returns:
# --------
# tempList : list
#   a list of the temperatures
#
proc ::xgs::util::optimalTempLadder {heatCapacity Tmin Tmax {prec 1}} {
    set N [optimalTempCount $heatCapacity $Tmin $Tmax]
    return [geometricTemperatures $Tmin $Tmax $N $prec]
}

# ::xgs::util::linearLambdas
#
# Return a lambda list with linear spacing:
#
# l_i = lmin + (lmax - lmin)(i-1)/(N-1)    i = 1, ..., N
#
# Arguments:
# ----------
# lmin : float
#   minimum lambda 
# lmax : float
#   maximum lambda 
# N : int
#   number of lambdas in ladder (including lmin and lmax)
# prec : int (optional, default: 2)
#   round lambda to "prec" decimal places
#
# Returns:
# --------
# lambdaList : list
#   a list of the lambda values
#
proc ::xgs::util::linearLambdas {lmin lmax N {prec 2}} {
    set lambdaList [list]
    set Scale [expr {pow(10, int($prec))}]
    set N [expr {int($N)}]
    set dl [expr {($lmax - $lmin) / ($N - 1)}]
    for {set i 1} {$i <= $N} {incr i} {
        set l [expr {$lmin + $dl*($i - 1)}]
        lappend lambdaList [expr {1.0*round($Scale*$l) / $Scale}]
    }
    return $lambdaList
}

# ::xgs::util::optimalLambdaCount
#
# Arguments:
# ----------
# d2fdl2 : float 
#   second derivative of the free energy wrt lambda
# lmin : float
#   minimum lambda 
# lmax : float
#   maximum lambda 
#
# Returns:
# --------
# Nopt : int
#   the number of lambda values
#
proc ::xgs::util::optimalLambdaCount {d2fdl2 lmin lmax} {
    return [expr {int(round(1 + 0.594*($lmax - $lmin)*sqrt(0.5*$d2fdl2)))}]
}

# ::xgs::util::optimalLambdaLadder 
#
# Arguments:
# ----------
# d2fdl2 : float 
#   second derivative of the free energy wrt lambda
# lmin : float
#   minimum lambda 
# lmax : float
#   maximum lambda 
# prec : int (optional, default: 2)
#   round lambda to "prec" decimal places
#
# Returns:
# --------
# lambdaList : list
#   a list of the lambda values
#
proc ::xgs::util::optimalLambdaLadder {d2fdl2 lmin lmax {prec 2}} {
    set N [optimalLambdaCount $d2fdl2 $lmin $lmax]
    return [linearLambdas $lmin $lmax $N $prec]
}

