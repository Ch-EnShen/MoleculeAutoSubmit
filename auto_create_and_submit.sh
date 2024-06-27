#!/bin/bash

# queue to number of cores mapping
declare -A queue_ncore
queue_ncore[mem256v3]=24 
queue_ncore[mem96]=32
queue_ncore[mem192]=32
queue_ncore[mem384]=32
queue_ncore[mem192v2]=32
queue_ncore[mem256v2]=32
queue_ncore[mem1T]=32

# declare some useful directories for later
mol_sub_dir=$(pwd)
dushin_path=$HOME/MOMAP_bak/momap/dushin
nacme_path=$HOME/MOMAP_bak/momap/nacme

####################### input PBS info #######################
# Read queue and set the corresponding ncore (default to mem192v2)
if [ ! -z "$2" ]; then
  queue="$2"
else
  read -p 'Queue? ' queue
fi

if [ -z $queue ]; then
  queue="mem192v2" #default
fi

ncore=${queue_ncore[$queue]}

# Reading of the nodes (default to 1)
if [ ! -z "$3" ]; then
  nnode="$3"
else
  read -p 'No. of Nodes? ' nnode
fi

if [ -z $nnode ]; then
  nnode="1" #default
fi

# Read the jobname (default to "mol_sub_default")
if [ ! -z "$4" ]; then
  jobname="$4"
else
  read -p 'Jobname? ' jobname
fi

if [ -z $jobname ]; then
  #jobname=$(cat qm.tpl | grep "sys_name" -A1 | tail -n 1 | cut -d'.' -f1)
  jobname="mol_sub_default"
fi

##############################################################

####################### define shell functions #######################
# modification of s0.gjf
function g16_sed_s0_opt(){ # $1: molecule_name, $2: functional
    x_temp=$queue
    x=${x_temp#mem}
    mem_new=${x%v*}
    sed -i "1c %chk=$1-$2-s0-opt-freq.chk" $1/$1-$2-s0-opt-freq.gjf
    sed -i "2c %mem=${mem_new}GB" $1/$1-$2-s0-opt-freq.gjf
    sed -i "3c %nprocshared=$ncore" $1/$1-$2-s0-opt-freq.gjf
}

# g16

function gen_g16wrapper_script(){
cat << EOF >> g16_wrapper.sh
#!/bin/bash

# Check if the correct number of arguments are provided
if [ \$# -lt 3 ]; then
    echo "Usage: \$0 inputfile queue jobname"
    exit 1
fi

# Read arguments
g16inputfile="\$1"
queue="\$2"  # Default to 'workq' if not provided
jobname="\$3"  # Default to input file name without extension

# Create an expect script to handle the interactive prompts
expect <<EOF
set timeout -1

# Start the g16sub script
spawn ~/bin/g16sub

# Handle the input file prompt
expect "G16 Input filename? " {
    send "\$g16inputfile\r"
}

# Handle the queue prompt
expect "Queue? " {
    send "\$queue\r"
}

# Handle the jobname prompt
expect "Jobname? " {
    send "\$jobname\r"
}

# Wait for the process to complete
expect eof
\EOF
EOF

chmod u+x g16_wrapper.sh
}

function error_termination(){ 
    # $1: mol_dir
    # $2: file_name.gjf
    
    # Normal termination or Error termination
    gjffile=$2
    filename="${gjffile%.gjf}" # convert .gjf to .log

    check=$(tail -n 10 ${1}/${filename}.log | grep 'Normal termination')
    if [ -z "$check" ]; then 
      echo "$2 process is Error terminated."
      return 1
    else
      echo "$2 process is Normal terminated."
    fi

    formchk ${filename}.chk
    mv $filename* $mol_dir/
    return 0
}

function wait_jobs(){
    for job in $(jobs -p)
    do
        wait $job
    done
    sleep 1
}

# dushin
function dushin_job(){ # for s0-s1 state
	  # $1: molecule_name, $2: functional
    mol_dir=${1}_${2}
    cp $mol_dir/$1-$2-*-opt-freq.fchk $dushin_path
    cp $mol_dir/$1-$2-*-opt-freq.log $dushin_path

    cd $dushin_path
    sed -i '1s/g16/g09/g' $dushin_path/$1-$2-s0-opt-freq.log
    sed -i '1s/g16/g09/g' $dushin_path/$1-$2-s1-opt-freq.log

    change8="2 1 'ground state' '$1-$2-s1-opt-freq.log'"
    change9="1 1 'excited state' '$1-$2-s0-opt-freq.log'"
    sed -i "8c $change8" $dushin_path/dushin.inp
    sed -i "9c $change9" $dushin_path/dushin.inp

    source $dushin_path/run 2> $1/dushin_res.txt

    cd $mol_sub_dir
    cp $dushin_path/dushin_cart.dus $mol_dir

    rm $dushin_path/$1*
    rm $dushin_path/dushin_cart.dus
}

# nacme
function nacme_job(){ 
    # only for s1 state
    # $1: mol_dir
    # $2: s1 file name

    cp $1/$2-opt-freq.log $nacme_path
    cp $1/$2-nacme.log $nacme_path

    cd $nacme_path
    sed -i "5s/get-nacme.*/get-nacme $2-nacme.log $2-opt-freq.log/g" $nacme_path/run
    source $nacme_path/run 2> $1/nacme_res.txt
    
    cd $mole_sub_dir/$2
    cp $nacme_path/nacme.out $1
    cp $nacme_path/nacme.log $1

    rm $nacme_path/$2*
    rm $nacme_path/nacme.out
    rm $nacme_path/nacme.log
}

function get_trans_dip_adiabatic_energy(){
    # $1: molecule_name, $2: functional
    python runscript_get_dipole_and_energy.py $1 $2
}

function cp_to_ICrate_workdir(){
    # $1: mol_dir
    account=${HOME#/home/}
    cp -r -T $1 /beegfs/hotpool/$account/Internal_Conversion/IC_Molecule/$1
}
###############################################################

####################### exection part #######################
functionals=(
    "b3lyp"
    "bmk"
    "pbe1pbe"
    "camb3lyp"
)

# wrapper for g16sub
gen_g16wrapper_script

for molecule_name in $(cat < molecule_name.txt)
do  
    for functional in "${functionals[@]}"; do
        # generate s1-opt-freq.gjf
        python generate_gjf_from_name.py $molecule_name $functional & 

        # remove undesired chars in molecule_name
        cleaned_molecule_name=$(echo "$molecule_name" | tr -d ',[]')
        mol_dir="${cleaned_molecule_name}_${functional}"

        # s0 opt freq
        g16_sed_s0_opt $cleaned_molecule_name $functional 
        s0_INPUTFILE="${mol_dir}/$cleaned_molecule_name-$functional-s0-opt-freq.gjf"
        ./g16_wrapper.sh $s0_INPUTFILE $queue $s0_INPUTFILE # g16sub s1-opt-freq.gjf
        wait_jobs
        if error_termination $mol_dir $s0_INPUTFILE; then # if Error termination, break the loop
            rm -r $mol_dir
            continue
        fi

        # s1 opt freq
        python runscript_mol_sub.py $cleaned_molecule_name opt-freq $functional & # generate s1-opt-freq.gjf
        s1_INPUTFILE="${mol_dir}/$cleaned_molecule_name-$functional-s1-opt-freq.gjf"
        ./g16_wrapper.sh $s1_INPUTFILE $queue $s1_INPUTFILE # g16sub s1-opt-freq.gjf
        wait_jobs
        if error_termination $mol_dir $s1_INPUTFILE; then # if Error termination, break the loop
            rm -r $mol_dir
            continue
        fi

        get_trans_dip_adiabatic_energy $cleaned_molecule_name $functional & # save as transition dipole, adiabatic energy as .npy file
        
        # nacme
        python runscript_mol_sub.py $cleaned_molecule_name nacme $functional & # generate s1-nacme.gjf
        nacme_INPUTFILE="${mol_dir}/$cleaned_molecule_name-$functional-s1-nacme.gjf"
        ./g16_wrapper.sh $nacme_INPUTFILE $queue $nacme_INPUTFILE # g16sub s1-nacme.gjf
        wait_jobs
        if error_termination $mol_dir $nacme_INPUTFILE; then # if Error termination, break the loop
            rm -r $mol_dir
            continue
        fi

        # dushin
        dushin_job $cleaned_molecule_name $functional & # generate dushin_cart.dus
        wait_jobs

        # final check: break the functional loop if Normal termination
        if !error_termination $mol_dir $s1_INPUTFILE; then # if Error termination, break the loop
            break
        fi
    done
done
###############################################################