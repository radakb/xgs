# numtcl.tcl
#
# Simple routines for operations on purely numerical lists. This is an 
# extremely weak attempt to regain some of my favorite functionality from numpy
# but with more Tcl-ish syntax.
#
# WARNING: USE AT YOUR OWN RISK, THERE IS GENERALLY VERY LITTLE ERROR CHECKING
# BEING PERFORMED. IF USED IMPROPERLY THESE ARE LIKELY TO FAIL WITH RATHER
# OPAQUE ERROR MESSAGES.
#

# tcl::mathfunc::lsum
#
# Sum the elements in a list.
#
proc tcl::mathfunc::lsum {list} {
    set sum 0.0
    foreach item $list {set sum [expr {$sum + $item}]}
    expr {$sum}
}

# tcl::mathfunc::lmean
#
# Average the elements in a list.
#
proc tcl::mathfunc::lmean {list} {
    set n [llength $list]
    set sum 0.0
    foreach item $list {set sum [expr {$sum + $item}]}
    expr {$sum/$n}
}

# lincr
#
# Increment a list (in place) by the given value. The modified list is also
# returned.
#
# - This supports the standard lindex syntax for nested lists (ala lindex).
# - NumPy style "broadcasting" is also supported, whereby a list plus a value
#   results in each element of the list being incremented and a list plus a
#   list results in element by element addition.
#   NB: This only works left to right (ala +=) and modifies the left list.
#
# Example:
#
# % set foo {0 2 1 7 3}
# % lincr foo 2
# 2 4 3 9 5
# % lincr foo {0 1 2 3 4}
# 2 5 5 12 9
#
# % set bar {0 1 {2 3 4} 5}
# % lincr bar 2 1
# 0 1 {3 4 5} 5
# % lincr bar 2 1 1
# 0 1 {3 5 5} 5
# 
proc lincr {list args} {
    upvar 1 $list List
    set value [lindex $args end]
    set index [lrange $args 0 end-1]
    
    set target [lindex $List {*}$index]
    set vlen [llength $value]
    set tlen [llength $target]

    if {$vlen == 1} {
        if {$tlen == 1} {
            # increment a single element by a single value
            lset List {*}$index [expr {$target + $value}]
        } elseif {$tlen > 1} {
            # increment a each element in a list by a single value
            set i 0 
            foreach item $target {
                if {[llength $item] != 1} {
                    error "lincr encountered a nested list during list + value\
                            operation"
                }
                lset List {*}$index $i [expr {$item + $value}]
                incr i
            }
        }
    } elseif {$vlen > 1} {
        if {$tlen != $vlen} {
            error "lincr cannot perform list + list operation with different\
                    lengths ($tlen != $vlen)"
        }
        set i 0
        foreach item1 $target item2 $value {
            lset List {*}$index $i [expr {$item1 + $item2}]
            incr i
        }
    }
    return $List
}

# lmultiply
#
# Multiple a list (in place) by the given value. The modified list is also
# returned.
#
# - This supports the standard lindex syntax for nested lists (ala lindex).
# - NumPy style "broadcasting" is also supported, whereby a list times a value
#   results in each element of the list being multiplied and a list times a
#   list results in element by element multiplication.
#   NB: This only works left to right (ala *=) and modifies the left list.
#
proc lmultiply {list args} {
    upvar 1 $list List
    set value [lindex $args end]
    set index [lrange $args 0 end-1]

    set target [lindex $List {*}$index]
    set vlen [llength $value]
    set tlen [llength $target]

    if {$vlen == 1} {
        if {$tlen == 1} {
            # increment a single element by a single value
            lset List {*}$index [expr {$target*$value}]
        } elseif {$tlen > 1} {
            # increment a each element in a list by a single value
            set i 0
            foreach item $target {
                if {[llength $item] != 1} {
                    error "lmultiply encountered a nested list during list *\
                            value operation"
                }
                lset List {*}$index $i [expr {$item*$value}]
                incr i
            }
        }
    } elseif {$vlen > 1} {
        if {$tlen != $vlen} {
            error "lmultiply cannot perform list * list operation with\
                    different lengths ($tlen != $vlen)"
        }
        set i 0
        foreach item1 $target item2 $value {
            lset List {*}$index $i [expr {$item1*$item2}]
            incr i
        }
    }
    return $List
}

# =============================================================================
# Probability, Statistics, and Random numbers
# =============================================================================
# tcl::mathfunc::uniform
#
# Generate a uniform random variable on the range [a, b).
#
# If only one argument is given, the range is [0, a) instead.
#
proc tcl::mathfunc::uniform {{a 1.0} {b ""}} {
    if {[string match $b ""]} {
        set b [expr {$a + 0.0}]
        set a 0.0
    }
    if {$a == $b} {
        return $a
    }
    if {$a > $b} {
        set tmp $a
        set a $b
        set b $tmp
    }
    expr {$a + ($b - $a)*rand()}
}

# tcl::mathfunc::randint
#
# Generate a random integer on the range [a, b].
#
proc tcl::mathfunc::randint {a b} {
    set a [expr {int($a)}]
    set b [expr {int($b)}]
    if {$a == $b} {
        return $a
    }
    if {$a > $b} {
        set tmp $a
        set a $b
        set b $tmp
    }
    expr {$a + int(($b - $a + 1)*uniform())}
}

# tcl::mathfunc::normal
#
# Generate a Gaussian (normal) random variable with the given mean and standard
# deviation (default to zero and one, respectively). 
#
# This uses the standard Box-Muller transformation of two uniform random 
# variables on (0, 1).
#
proc tcl::mathfunc::normal {{mu 0.0} {sigma 1.0}} {
    set two_pi [expr {2*3.141592653589793}]
    set u1 [expr {uniform()}]
    set u2 [expr {uniform()}]
    # The transformation produces two random variates - either of these can be
    # used (or both).
    expr {$mu + $sigma*sqrt(-2*log($u1))*cos($two_pi*$u2)}
    # expr {$mu + $sigma*sqrt(-2*log($u1))*sin($two_pi*$u2)}
}

# choice
#
# Choose a random element from a list. The index of the element is also
# returned.
#
# The default is to uniformly weight the elements, but (unnormalized)
# non-uniform weights may also be provided.
#
# Arguments:
# ----------
# list : list
#   list to choose an element from
# weights : list of floats (optional)
#   The probability weight of each choice. These need not be normalized.
#
# Returns:
# --------
# element
#   the selected list element
# index : int
#   the index of the selected element
#
proc choice {list {weights {}}} {
    # Uniform weights is really easy.
    if {![llength $weights]} {
        set j [expr {int(rand()*[llength $list])}]
        return [list [lindex $list $j] $j]
    }
    # Non-uniform weighting is more involved.
    if {[llength $list] != [llength $weights]} {
        error "Mismatch in list/weights in choice!"
    }
    set WeightSum [expr {lsum($weights)}]
    set Rand [expr {$WeightSum*rand()}]
    set j 0
    foreach Weight $weights {
        set Rand [expr {$Rand - $Weight}]
        if {$Rand <= 0.0} {
            return [list [lindex $list $j] $j]
        }
        incr j
    }
    # This should never be reached.
    error "Something bad happened in choice!"
}

# linearFit
#
# Perform a (possibly weighted) linear least squares fit.
#
# Arguments:
# ----------
# x : list
#   the x-values to fit (must be same size as y)
# y : list
#   the y-values to fit (must be same size as x)
# yerr : list (optional)
#   uncertainty in y, will be used as inverse weights if specified
# intrcpt : float (optional)
#   the fixed y-intercept to use during fitting
#
# Returns:
# --------
# slp : float
#   the fitted slope
# intrcpt : float
#   the fitted (or fixed) intercept
# slp_err : float
#   the uncertainty in the slope
# intrcpt_err : float
#   the uncertainty in the intercept (non-zero even for fixed slope)
# R2 : float
#   the squared coefficient of regression
#
proc linearFit {x y {yerr {}} {intrcpt {}}} {
    if {[llength $x] != [llength $y]} {
        error "The x and y lists are not of the same length!"
    }
    set n [llength $x]
    set dof [expr {$n - 2}]
    if {[llength $yerr] == 0} {
        set w [lrepeat $n 1.0]
    } else {
        set w [list]
        foreach yerri $yerr {
            lappend w [expr {1.0/$yerri}]
        }
    }
    if {[llength $w] != $n} {
        error "The x/y and weight lists are not of the same length!"
    }
    set fixed_intrcpt [expr {$intrcpt == "" ? 0 : 1}]
    if {$fixed_intrcpt} {
        lincr y [lrepeat $n [expr {-1*$intrcpt}]]
        incr dof 1
    }
    # Normalize the data by the sample weights.
    lmultiply x $w
    lmultiply y $w
    # Compute the covariance matrix elements.
    set xmean [expr {lmean($x)}]
    set ymean [expr {lmean($y)}]
    set dx $x
    set dy $y
    lincr dx [lrepeat $n [expr {-1*$xmean}]]
    lincr dy [lrepeat $n [expr {-1*$ymean}]]
    set xvar 0.0
    set yvar 0.0
    set xyvar 0.0
    foreach dxi $dx dyi $dy {
        set xvar [expr {$xvar + $dxi**2}]
        set yvar [expr {$yvar + $dyi**2}]
        set xyvar [expr {$xyvar + $dxi*$dyi}]
    }
    set xvar [expr {$xvar / $n}]
    set yvar [expr {$yvar / $n}]
    set xyvar [expr {$xyvar / $n}]
    # Shift the data if the intercept is not being calculated.
    if {$fixed_intrcpt} {
        set xvar [expr {$xvar + $xmean**2}]
        set yvar [expr {$yvar + $ymean**2}]
        set xyvar [expr {$xyvar + $xmean*$ymean}]
    }
    set slp [expr {$xyvar / $xvar}]
    # Compute the intercept if not already specified.
    if {!$fixed_intrcpt} {
        set intrcpt [expr {$ymean - $slp*$xmean}]
    }
    # Compute the correlation coefficient and parameter variances.
    set R2 [expr {$xyvar**2 / ($xvar*$yvar)}]
    set res_var 0.0
    foreach xi $x yi $y {
        set res_var [expr {$res_var + ($slp*$xi + $intrcpt - $yi)**2}]
    }
    set res_var [expr {$res_var / $dof}]
    lmultiply dx $dx
    set slp_var [expr {$res_var / lsum($dx)}] ;# actually divide by dx**2
    lmultiply x $x
    set intrcpt_var [expr {$slp_var*lmean($x)}] ;# actually divide by x**2
    set slp_err [expr {sqrt($slp_var/$dof)}]
    set intrcpt_err [expr {sqrt($intrcpt_var/$dof)}]
    return [list $slp $intrcpt $slp_err $intrcpt_err $R2]
}

