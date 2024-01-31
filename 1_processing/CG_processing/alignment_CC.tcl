set sellist {}

set nChains 10
set atomsPerMol 59
set fin [expr $atomsPerMol-1]
mol new vs.lammpstrj waitfor all autobonds off
set nf [molinfo top get numframes]
#set nf 3
set fini [expr $fin - 3]
set AAref [atomselect top "index 0 to $fin" frame 0]
set com [measure center $AAref weight mass]
$AAref moveby [vecscale -1.0 $com]
set reference [atomselect top "index 9 to $fini" frame 0]
for {set i 0} {$i<$nChains} {incr i} {
	puts "starting chain $i"

	for {set f 0} {$f<$nf} {incr f} {
		#puts $f
		set rangei [expr $i*$atomsPerMol+9]
		set rangef [expr $i*$atomsPerMol+$atomsPerMol-1-3]
		
		set init [expr $i*$atomsPerMol]
		set end [expr $i*$atomsPerMol+$atomsPerMol-1]
		set target [atomselect top "index $rangei to $rangef" frame $f]
		set targetall [atomselect top "index $init to $end" frame $f]
		set M [measure fit $target $reference]
		$targetall move $M
		#lappend sellist $complete_selection
		#$reference delete
		$target delete
		$targetall delete
		#$targetall_2 delete
	}
}
#set molall [::TopoTools::selections2mol $sellist]
 
for {set i 0} {$i<$nChains} {incr i} {
	set init [expr $i*$atomsPerMol]
	set end [expr $i*$atomsPerMol+$atomsPerMol-1]
	set tl [atomselect top "index $init to $end"]
	animate write lammpstrj CA_monomers${i}.lammpstrj sel $tl
}
#animate write lammpstrj vs_m.lammpstrj $sellist
exit
