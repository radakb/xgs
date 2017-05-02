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
# Increment an element of a list by the given value.
# Indexing uses the standard syntax for nested lists, similar to lindex.
#
proc lincr {list args} {
    upvar 1 $list List
    set value [lindex $args end]
    set index [lrange $args 0 end-1]
    lset List {*}$index [expr {[lindex $List {*}$index] + $value}]
    return
}

# llincr
#
# Increment each element of list1 by the elements in list2.
#
proc llincr {list1 args} {
    upvar 1 $list1 List1
    set list2 [lindex $args end]
    set index [lrange $args 0 end-1]
    set i 0
    foreach item1 [lindex $List1 {*}$index] item2 $list2 {
        lset List1 {*}$index $i [expr {$item1 + $item2}]
        incr i
    }
    return
}

# lmultiply
#
# Multiple each element of a list by the given value.
#
proc lmultiply {list args} {
    upvar 1 $list List
    set value [lindex $args end]
    set index [lrange $args 0 end-1]
    set i 0
    foreach item [lindex $List {*}$index] {
        lset List {*}$index $i [expr {$item*$value}]
        incr i
    }
    return
}

# ldot
#
# Multiply a list (in place) by another list (i.e. take the dot product).
#
proc ldot {list1 list2} {
    upvar 1 $list1 List1
    set i 0
    foreach item1 $List1 item2 $list2 {
        lset List1 $i [expr {$item1*$item2}]
        incr i
    }
    return
}

# =============================================================================
# Probability, Statistics, and Random numbers
# =============================================================================
# Generate a Gaussian (normal) random variable with the given mean and standard
# deviation (default to zero and one, respectively). 
#
# This uses the standard Box-Muller transformation of two uniform random 
# variables on (0, 1).
#
proc tcl::mathfunc::normal {{mu 0.0} {sigma 1.0}} {
    set two_pi [expr {2*3.141592653589793}]
    set u1 [expr {rand()}]
    set u2 [expr {rand()}]
    # The transformation produces two random variates - either of these can be
    # used (or both).
    expr {$mu + $sigma*sqrt(-2*log($u1))*cos($two_pi*$u2)}
    # expr {$mu + $sigma*sqrt(-2*log($u1))*sin($two_pi*$u2)}
}

# choice
#
# Choose a random element from a list with equal probability.
#
# Arguments:
# ----------
# list : list
#    list to choose an element from
#
proc choice {list} {
    return [lindex $list [expr {int(rand()*[llength $list])}]]
}

# weightedChoice
#
# Given a set of probability weights, return the index of a given outcome with 
# the correct probability.
#
# Arguments:
# ----------
# weights : list of floats
#   The probability weight of each choice. These need not be normalized.
# weightSum : float (optional, calculated if not specified)
#   The sum of all entries in weights. This will not function correctly if this
#   is incorrectly specified!
#
proc weightedChoice {weights {weightSum {}}} {
    if {[llength $weightSum]} {
        set WeightSum $weightSum
    } else {
        set WeightSum [expr {lsum($weights)}]
    }
    set Rand [expr {$WeightSum*rand()}]
    set j 0
    foreach Weight $weights {
        set Rand [expr {$Rand - $Weight}]
        if {$Rand <= 0.0} {
            return $j
        }
        incr j
    }
    # This should never be reached.
    error "Something bad happened in weightedChoice!"
}

# linearFit
#
# Perform a (weighted) linear least squares fit.
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
        llincr y [lrepeat $n [expr {-1*$intrcpt}]]
        incr dof 1
    }
    # Normalize the data by the sample weights.
    ldot x $w
    ldot y $w
    # Compute the covariance matrix elements.
    set xmean [expr {lmean($x)}]
    set ymean [expr {lmean($y)}]
    set dx $x
    set dy $y
    llincr dx [lrepeat $n [expr {-1*$xmean}]]
    llincr dy [lrepeat $n [expr {-1*$ymean}]]
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
    ldot dx $dx
    set slp_var [expr {$res_var / lsum($dx)}] ;# actually divide by dx**2
    ldot x $x
    set intrcpt_var [expr {$slp_var*lmean($x)}] ;# actually divide by x**2
    set slp_err [expr {sqrt($slp_var/$dof)}]
    set intrcpt_err [expr {sqrt($intrcpt_var/$dof)}]
    return [list $slp $intrcpt $slp_err $intrcpt_err $R2]
}
