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
  read -p 'Queue? (default to mem192v2)' queue
fi

if [ -z $queue ]; then
  queue="mem192v2" #default
fi

ncore=${queue_ncore[$queue]}
##############################################################

####################### define shell functions #######################
# start job
function start_job(){
    # $1: cleaned molecule name
    # $2: jobname
cat << EOF > mol_calc_$1.sh
#!/bin/bash
#PBS -N $2
#PBS -o $2.out
#PBS -e $2.err
#PBS -q $queue
#PBS -l nodes=1:ppn=$ncore

cd \$PBS_O_WORKDIR/
module purge
module load g16-a03
module load anaconda3.8

start=\$(date "+%s.%N")

pip install pubchempy

EOF
}

# modification of s0.gjf
function g16_sed_s0_opt(){ 
    # $1: molecule_name
    x_temp=$queue
    x=${x_temp#mem}
    mem_new=${x%v*}
    sed -i "1c %chk=${1}-\${functional}-s0-opt-freq.chk" ${1}_\${functional}/${1}-\${functional}-s0-opt-freq.gjf
    sed -i "2c %mem=${mem_new}GB" ${1}_\${functional}/$1-\${functional}-s0-opt-freq.gjf
    sed -i "3c %nprocshared=$ncore" ${1}_\${functional}/$1-\${functional}-s0-opt-freq.gjf
}

function python_create_s0_gjf(){
    # $1: molecule name
    # $2: cleaned molecule name
cat << EOF >> mol_calc_$2.sh
    python generate_gjf_from_name.py $1 \$functional

    for job in \$(jobs -p)
    do
        wait \$job
    done

EOF
}

# g16
function g16job(){
    # $1: clean_molecule_name
    # $2: state (s0, s1)
    # $3: job_type (opt-freq, nacme)
cat << EOF >> mol_calc_$1.sh
    INPUTFILE=${1}_\${functional}/$1-\${functional}-$2-$3.gjf

    # Setup G09 runtime environment
    OUTPUTFILE=\${INPUTFILE%.*}.log
    export GAUSS_SCRDIR=/tmp/\$USER/\$RANDOM
    mkdir -p \$GAUSS_SCRDIR

    echo "======================================================="
    echo "Working directory is" \$PBS_O_WORKDIR
    echo "Starting on \$(hostname) at \$(date)"


    # Run g16 job
    cd \$PBS_O_WORKDIR/
    time g16 < \$INPUTFILE > \$OUTPUTFILE

    # Clean scratch file once finished
    rm -rf \$GAUSS_SCRDIR

    for job in \$(jobs -p)
    do
        wait \$job
    done
EOF
}

function python_mol_sub(){
    # $1: clean_molecule_name
    # $2: job_type (opt-freq, nacme)
cat << EOF >> mol_calc_$1.sh
    python runscript_mol_sub.py $1 $2 \$functional

    for job in \$(jobs -p)
    do
        wait \$job
    done
EOF
}

# dushin
function dushin_job(){ # for s0-s1 state
  # $1: cleaned molecule name
cat << EOF >> mol_calc_$1.sh
    cp \$PBS_O_WORKDIR/${1}_\${functional}/$1-\${functional}-*-opt-freq.fchk $dushin_path
    cp \$PBS_O_WORKDIR/${1}_\${functional}/$1-\${functional}-*-opt-freq.log $dushin_path

    cd $dushin_path
    sed -i '1s/g16/g09/g' $dushin_path/$1-\${functional}-s0-opt-freq.log
    sed -i '1s/g16/g09/g' $dushin_path/$1-\${functional}-s1-opt-freq.log

    change8="2 1 'ground state' '$1-\${functional}-s1-opt-freq.log'"
    change9="1 1 'excited state' '$1-\${functional}-s0-opt-freq.log'"
    sed -i "8c \$change8" $dushin_path/dushin.inp
    sed -i "9c \$change9" $dushin_path/dushin.inp

    source $dushin_path/run 2> \$PBS_O_WORKDIR/${1}_\${functional}/dushin_res.txt

    cp $dushin_path/dushin_cart.dus \$PBS_O_WORKDIR/${1}_\${functional}/

    rm $dushin_path/$1*
    rm dushin_cart.dus
EOF
}

# nacme
function nacme_job(){ # only for s1 state
    # $1: cleaned molecule name
cat << EOF >> mol_calc_$1.sh
    cp \$PBS_O_WORKDIR/${1}_\${functional}/$1-\${functional}-s1-opt-freq.log $nacme_path
    cp \$PBS_O_WORKDIR/${1}_\${functional}/$1-\${functional}-s1-nacme.log $nacme_path

    cd $nacme_path
    sed -i "5s/get-nacme.*/get-nacme $1-\${functional}-s1-nacme.log $1-\${functional}-s1-opt-freq.log/g" $nacme_path/run

    source $nacme_path/run 2> \$PBS_O_WORKDIR/${1}_\${functional}/nacme_res.txt
    cp $nacme_path/nacme.out \$PBS_O_WORKDIR/${1}_\${functional}/
    cp $nacme_path/nacme.log \$PBS_O_WORKDIR/${1}_\${functional}/

    rm $nacme_path/$1*
    rm nacme.out
    rm nacme.log
    
    for job in \$(jobs -p)
    do
        wait \$job
    done
    
    cd \$PBS_O_WORKDIR
EOF
}

function get_trans_dip_adiabatic_energy(){
    # $1: cleaned molecule name
cat << EOF >> mol_calc_$1.sh
    python runscript_get_dipole_and_energy.py $1 \${functional}
    for job in \$(jobs -p)
    do
        wait \$job
    done

    sleep 1
EOF
}

function cp_to_ICrate_workdir(){
    # $1: cleaned molecule name
cat << EOF >> mol_calc_$1.sh
    account=\${HOME#/home/}
    cd \$PBS_O_WORKDIR
    cp -r -T ${1}_\${functional} /beegfs/hotpool/\$account/Internal_Conversion/IC_Molecule/${1}_\${functional}
EOF
}
###############################################################

#### check part
function normal_termination(){
    # Normal termination or Error termination

    # $1: clean_molecule_name
    # $2: state (s0, s1)
    # $3: job_type (opt-freq, nacme)
cat << EOF >> mol_calc_$1.sh
    # mol_dir="${1}_\${functional}"
    # filename=${1}-\${functional}-${2}-${3}
    
    formchk ${1}-\${functional}-${2}-${3}.chk
    mv ${1}-\${functional}-${2}-${3}* ${1}_\${functional}/
    
    if [ "$2" == "s0" ]; then
      check=\$(tail -n 1 "${1}_\${functional}/${1}-\${functional}-${2}-${3}.log" | grep 'Normal termination')
      if [ -z "\$check" ]; then 
        echo "${1}-\${functional}-${2}-${3} process is Error terminated."
        continue
      else
        echo "${1}-\${functional}-${2}-${3} process is Normal terminated."
      fi
    elif [ "$2" == "s1" ] && [ "$3" == "nacme" ]; then
      check=\$(tail -n 1 "${1}_\${functional}/${1}-\${functional}-${2}-${3}.log" | grep 'Normal termination')
      if [ -z "\$check" ]; then
        echo "${1}-\${functional}-${2}-${3} process is Error terminated."
      else
        echo "${1}-\${functional}-${2}-${3} process is Normal terminated."
        break
      fi
    elif [ "$2" == "s1" ]; then
      check=\$(tail -n 1 "${1}_\${functional}/${1}-\${functional}-${2}-${3}.log" | grep 'Normal termination')
      if [ -z "\$check" ]; then
        echo "${1}-\${functional}-${2}-${3} process is Error terminated."
        continue
      else
        echo "${1}-\${functional}-${2}-${3} process is Normal terminated."
      fi
    else
      echo "Unknown state: $3"
      continue
    fi
EOF
}

function end_script(){
cat << EOF >> rate_caln.sh
end=\$(date "+%s.%N")

echo "Job Ended at \$(date)"
echo '======================================================='
#mpdallexit
EOF
}
####

####################### execution part #######################
functionals=(
    "b3lyp"
    "bmk"
    "pbe1pbe"
    "camb3lyp"
)

for molecule_name in $(cat < molecule_name.txt)
do
    cleaned_molecule_name=$(echo "$molecule_name" | tr -d ',[]')
    
    # setup for building molecular directory
    for functional in "${functionals[@]}"
    do
      python generate_gjf_from_name.py $molecule_name $functional
    done
    
    # setup for qsub
    start_job $cleaned_molecule_name $cleaned_molecule_name
cat << EOF >> mol_calc_$cleaned_molecule_name.sh
functionals=(
    "b3lyp"
    "bmk"
    "pbe1pbe"
    "camb3lyp"
)

for functional in "\${functionals[@]}"
do
EOF
    # python_create_s0_gjf $molecule_name $cleaned_molecule_name # generate s0-opt-freq.gjf
    g16job $cleaned_molecule_name s0 opt-freq # calculate s0-opt-freq
    normal_termination $cleaned_molecule_name s0 opt-freq # check whether process error for s0-opt-freq
    python_mol_sub $cleaned_molecule_name opt-freq # generate s1-opt-freq.gjf
    g16job $cleaned_molecule_name s1 opt-freq # calculate s1-opt-freq
    normal_termination $cleaned_molecule_name s1 opt-freq # check whether process error for s1-opt
    get_trans_dip_adiabatic_energy $cleaned_molecule_name # save as transition dipole, adiabatic energy as .npy file
    python_mol_sub $cleaned_molecule_name nacme # generate s1-nacme.gjf
    g16job $cleaned_molecule_name s1 nacme # calculate s1-nacme
    dushin_job $cleaned_molecule_name # generate dushin_cart.dus
    nacme_job $cleaned_molecule_name # generate nacme.out
    normal_termination $cleaned_molecule_name s1 nacme # check whether process error for s1-nacme

cat << EOF >> mol_calc_$cleaned_molecule_name.sh
done
EOF

end_script

qsub mol_calc_$cleaned_molecule_name.sh
rm mol_calc_$cleaned_molecule_name.sh

done
###############################################################