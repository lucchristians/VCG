
mol new imflex.lammpstrj waitfor all autobonds off
set nf [molinfo top get numframes]
set AAref [atomselect top "all" frame 0]
set atomsPerMol [$AAref num]
set fin [expr $atomsPerMol-1]
set fini [expr $fin - 3]
set com [measure center $AAref weight mass]
$AAref moveby [vecscale -1.0 $com]
set reference [atomselect top "index 5 to $fini" frame 0]

for {set f 1} {$f<$nf} {incr f} {
	set rangei [expr 5]
	set rangef [expr $atomsPerMol-1-3]
	set init [expr 0]
	set end [expr $atomsPerMol-1]
	set target [atomselect top "index $rangei to $rangef" frame $f]
	set targetall [atomselect top "all" frame $f]
	set M [measure fit $target $reference]
	$targetall move $M
	$target delete
	$targetall delete
}

set tl [atomselect top "all"]
animate write lammpstrj imflex_a.lammpstrj sel $tl
exit
