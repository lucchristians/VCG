set sellist {}

set nChains 5
set atomsPerMol 54
set fin [expr $atomsPerMol-1]
mol new AA_novs.lammpstrj waitfor all autobonds off
set nf [molinfo top get numframes]
#set nf 3
set fini [expr $fin - 3]
set AAref [atomselect top "index 0 to $fin 540" frame 0]
set com [measure center $AAref weight mass]
$AAref moveby [vecscale -1.0 $com]
set reference [atomselect top "index 9 to $fini" frame 0]
for {set i 0} {$i<$nChains} {incr i} {
	puts "starting chain $i"
	set vs1_index [expr $i * 2 + 540]
	set vs2_index [expr $i * 2 + 541]
	puts $vs1_index
	puts $vs2_index

	for {set f 0} {$f<$nf} {incr f} {
		#puts $f
		set rangei [expr $i*$atomsPerMol+9]
		set rangef [expr $i*$atomsPerMol+$atomsPerMol-1-3]
		set rangei_2 [expr 49+$i*$atomsPerMol]
		set rangef_2 [expr 53+$i*$atomsPerMol]
		set targeti_2 [expr 319+$i*$atomsPerMol]
		set targetf_2 [expr 323+$i*$atomsPerMol]
		
		set init [expr $i*$atomsPerMol]
		set end [expr $i*$atomsPerMol+$atomsPerMol-1]
		set target [atomselect top "index $rangei to $rangef" frame $f]
		set targetall [atomselect top "index $init to $end $vs1_index" frame $f]
		set M [measure fit $target $reference]
		$targetall move $M
		$target num
		set target2 [atomselect top "index $targeti_2 to $targetf_2" frame $f]
		#puts [$target2 num]
		#set targetall_2 [atomselect 0 "index $targeti_2 to $targetf_2 $vs2_index" frame $f]
		
		set reference2 [atomselect top "index $rangei_2 to $rangef_2" frame $f]
		set vs2 [atomselect top "index ${vs2_index}" frame $f]
		#puts [$reference2 num]
		set M2 [measure fit $target2 $reference2]
		#puts $M2
		$vs2 move $M2 
		#set complete_selection [atomselect 0 "$init to $end $vs1_index $vs2_index"]
		#lappend sellist $complete_selection
		#$reference delete
		$target delete
		$target2 delete
		$reference2 delete
		$targetall delete
		$vs2 delete
		#$targetall_2 delete
	}
}
#set molall [::TopoTools::selections2mol $sellist]
 
for {set i 0} {$i<$nChains} {incr i} {
	set init [expr $i*$atomsPerMol]
	set end [expr $i*$atomsPerMol+$atomsPerMol-1]
	set vs1_index [expr $i * 2 + 540]
	set vs2_index [expr $i * 2 + 541]
	set tl [atomselect top "index $init to $end $vs1_index $vs2_index"]
	animate write lammpstrj CA_m${i}.lammpstrj sel $tl
}
#animate write lammpstrj vs_m.lammpstrj $sellist
exit
