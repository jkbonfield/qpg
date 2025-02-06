#!/usr/bin/wish

#wm withdraw .

proc query_len {cigar} {
    set qlen 0
    foreach {cig len op} [regexp -all -inline {([0-9]+)(.)} $cigar] {
	switch $op {
	    M -
	    S -
	    H -
	    I {incr qlen $len}
	}
    }
    return $qlen
}

set width  1000
set height 1000
set c [canvas .c -width $width -height $height]
pack $c -fill both -expand 1

set colours [list blue yellow red2 green pink4 purple]
set cidx 0

set fd [open [lindex $argv 0]]
set max_qpos 1
set max_rpos 1

# For multiple contigs, we concatenate to a single position
set total_qlen 0
set last_contig ""
set last_qlen 0

while {[gets $fd line] != -1} {
    set colour [lindex $colours $cidx]
    if {[incr cidx] == [llength $colours]} {
	set cidx 0
    }

    if {[string match "@*" $line]} continue
    set f [split $line "\t"]

    if {$last_contig != "" && [lindex $line 0] != $last_contig} {
	puts "New contig [lindex $line 0], last_len $last_qlen"
	incr total_qlen $last_qlen
	$c create line -1e10 $total_qlen 1e10 $total_qlen
    }
    set last_contig [lindex $line 0]

    set cigar [lindex $line 5]
    set rpos [lindex $line 3]
    set rev [expr {[lindex $line 1] & 16}]
    #set qpos [expr {$rev ? [string length [lindex $line 9]] : 0}]
    set last_qlen [query_len $cigar]
    set qpos [expr {$rev ? $last_qlen : 0}]
    incr qpos $total_qlen
    #incr qpos [lindex $line 3]; # POS
    puts -nonewline "[lindex $line 0] $rev    $rpos $qpos // "
    set segment 0
    set x $rpos
    set y $qpos

    if {$max_rpos < $rpos} {
	set max_rpos $rpos
    }
    if {$max_qpos < $qpos} {
	set max_qpos $qpos
    }

    foreach {cig len op} [regexp -all -inline {([0-9]+)(.)} $cigar] {
	set r1 $rpos
	set q1 $qpos
	set qlen [expr {$rev ? -$len : $len}];
	switch $op {
	    X -
	    = -
	    M {
		set rpos [expr {$rpos+$len}]
		set qpos [expr {$qpos+$qlen}]
	    }
	    D {
		set rpos [expr {$rpos+$len}]
	    }
	    S -
	    H -
	    I {
		set qpos [expr {$qpos+$qlen}]
	    }
	}

	#set c2 [tk::Darken $colour [expr {80+($segment%2)*40}]]
	if {$op != "S" && $op != "H"} {
	    #puts "\t$len $op\t$r1 $q1 ... $rpos $qpos"
	    $c create line $r1 $q1 $rpos $qpos \
		-width 5 \
		-fill $colour
	}
    }
    puts "$rpos $qpos // $colour"

    if {$max_rpos < $rpos} {
	set max_rpos $rpos
    }
    if {$max_qpos < $qpos} {
	set max_qpos $qpos
    }
}

.c scale all 0 0 \
    [expr {double($width)/$max_rpos}] \
    [expr {double($height)/$max_qpos}]
.c scale all 600 600 0.9 0.9
#exit


