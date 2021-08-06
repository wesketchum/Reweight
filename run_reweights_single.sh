#!/bin/sh

for folder in /Data/LDMX/LDMX-eN/weights_28Jul2021/*; do
    echo $folder
    input_file=$(ls $folder/gntp.0.ghep.root)
    echo $input_file
    output_loc=$folder

    #input_file=/Data/LDMX/gntp.0.ghep.root
    #output_loc=/Data/LDMX
    seed=1000
    ntwk=5
    min_twk=-2
    max_twk=2

    knobs=(
	"FrAbs_N"
	"FrInel_N"
	"FrPiProd_N"
	"FrCEx_N"
	"FrAbs_pi"
	"FrInel_pi"
	"FrPiProd_pi"
	"FrCEx_pi"
)
#    "MFP_N"
#    "MFP_pi"
#    "FormZone"
#)
#    "AGKYxF1pi"
#    "AGKYpT1pi"
#)

#grwght1p -f /Data/LDMX/gntp.0.ghep.root -s FrAbs_N -t 5 --min-tweak -2 --max-tweak 2 --seed 1000 -p 11 -o /Data/LDMX/weights_FrAbs_N_all.root >& weights_FrAbs_N_all.log &

    for k in "${knobs[@]}"; do
	echo "Launch knob $k"
	out_filename="${output_loc}/weights_${k}_FSIFix.root"
	log_file="${output_loc}/weights_${k}.log"
	echo $out_filename
	grwght1p -f $input_file -s $k -t $ntwk --min-tweak $min_twk --max-tweak $max_twk --seed $seed -p 11 -o $out_filename >& $log_file &
	sleep 1
    done

    echo "Done with $folder"


    #exit 0

    sleep 3600



done
