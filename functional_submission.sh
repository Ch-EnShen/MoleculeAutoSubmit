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

####################### input PBS info #######################
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

# Read the jobname (default to "IC_default")
if [ ! -z "$4" ]; then
  jobname="$4"
else
  read -p 'Jobname? ' jobname
fi

if [ -z $jobname ]; then
  #jobname=$(cat qm.tpl | grep "sys_name" -A1 | tail -n 1 | cut -d'.' -f1)
  jobname="mol_sub_default"
fi

# Read the functional (default to "b3lyp")
if [ ! -z "$5" ]; then
  functional="$5"
else
  read -p 'Functional? ' functional
fi

if [ -z $functional ]; then
  functional="b3lyp"
fi
##############################################################

####################### create file s0.gjf #######################
function python_create_s0_gjf(){
cat << EOF >> mol_calc.sh
    python generate_gjf_from_name.py \$molecule_name $functional &

    for job in \$(jobs -p)
    do
        wait \$job
    done

EOF
}

####################### modification of s0.gjf #######################
function g16_sed_s0_opt(){
    # $1: molecule_name
    x_temp=$queue
    x=${x_temp#mem}
    mem_new=${x%v*}
    sed -i "1c %chk=${1}-${functional}-s0-opt-freq.chk" $1/$1-$functional-s0-opt-freq.gjf
    sed -i "2c %mem=${mem_new}GB" $1/$1-$functional-s0-opt-freq.gjf
    sed -i "3c %nprocshared=$ncore" $1/$1-$functional-s0-opt-freq.gjf
}
######################################################################

##################### g16, dushin, nacme ##################### 
dushin_path=$HOME/MOMAP_bak/momap/dushin
nacme_path=$HOME/MOMAP_bak/momap/nacme

function g16job(){
    # $1: state (s0, s1)
    # $2: job_type (opt-freq, nacme)
cat << EOF >> mol_calc.sh
    INPUTFILE=\$molecule_name/\$molecule_name-$functional-$1-$2.gjf

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

    # Check whether process error
    check=\$(tail -n 1 \$molecule_name/\$molecule_name-$functional-$1-$2.log | grep 'Normal termination')
    if [ -z "\$check" ]; then 
      echo "$1-$2 process is not normally terminated."
      continue
    fi

    formchk \$molecule_name-$functional-$1-$2.chk
    mv \$molecule_name-$functional-$1-$2* \$molecule_name/

EOF
}

function python_mol_sub(){
    # $1: job_type
cat << EOF >> mol_calc.sh
    python runscript_mol_sub.py \$molecule_name $1 $functional&

    for job in \$(jobs -p)
    do
        wait \$job
    done
EOF
}

function dushin_job(){ # for s0-s1 state
	# $1: functional
cat << EOF >> mol_calc.sh
    cp \$PBS_O_WORKDIR/\$molecule_name/\$molecule_name-$functional-*-opt-freq.fchk $dushin_path
    cp \$PBS_O_WORKDIR/\$molecule_name/\$molecule_name-$functional-*-opt-freq.log $dushin_path

    cd $dushin_path
    sed -i '1s/g16/g09/g' $dushin_path/\$molecule_name-$functional-s0-opt-freq.log
    sed -i '1s/g16/g09/g' $dushin_path/\$molecule_name-$functional-s1-opt-freq.log

    change8="2 1 'ground state' '\$molecule_name-$functional-s1-opt-freq.log'"
    change9="1 1 'excited state' '\$molecule_name-$functional-s0-opt-freq.log'"
    sed -i "8c \$change8" $dushin_path/dushin.inp
    sed -i "9c \$change9" $dushin_path/dushin.inp

    source $dushin_path/run 2> \$PBS_O_WORKDIR/\$molecule_name/dushin_res.txt

    cp $dushin_path/dushin_cart.dus \$PBS_O_WORKDIR/\$molecule_name/

    rm $dushin_path/\$molecule_name*
    rm dushin_cart.dus
EOF
}

function nacme_job(){ # only for s1 state
cat << EOF >> mol_calc.sh
    cp \$PBS_O_WORKDIR/\$molecule_name/\$molecule_name-$functional-s1-opt-freq.log $nacme_path
    cp \$PBS_O_WORKDIR/\$molecule_name/\$molecule_name-$functional-s1-nacme.log $nacme_path

    cd $nacme_path
    sed -i "5s/get-nacme.*/get-nacme \$molecule_name-$functional-s1-nacme.log \$molecule_name-$functional-s1-opt-freq.log/g" $nacme_path/run

    source $nacme_path/run 2> \$PBS_O_WORKDIR/\$molecule_name/nacme_res.txt
    cp $nacme_path/nacme.out \$PBS_O_WORKDIR/\$molecule_name/
    cp $nacme_path/nacme.log \$PBS_O_WORKDIR/\$molecule_name/

    rm $nacme_path/\$molecule_name*
    rm nacme.out
    rm nacme.log
    
    for job in \$(jobs -p)
    do
        wait \$job
    done

    end=\$(date "+%s.%N")
    runtime=\$(echo "\$end - \$start" | bc -lq)

    echo "Job Ended at \$(date)"
    echo "Estimated Runtime: \$runtime seconds."
    echo '======================================================='
    #mpdallexit
EOF
}

function get_trans_dip_adiabatic_energy(){
cat << EOF >> mol_calc.sh
    python runscript_get_dipole_and_energy.py \$molecule_name $functional
    for job in \$(jobs -p)
    do
        wait \$job
    done

    sleep 1
EOF
}

function cp_to_ICrate_workdir(){
cat << EOF >> mol_calc.sh
    account=\${HOME#/home/}
    cd \$PBS_O_WORKDIR
    cp -r -T \$molecule_name /beegfs/hotpool/\$account/Internal_Conversion/IC_Molecule/\$molecule_name_$functional
EOF
}
################################################################

function start_job(){
cat << EOF > mol_calc.sh
#!/bin/bash
#PBS -N $jobname
#PBS -o $jobname.out
#PBS -e $jobname.err
#PBS -q $queue
#PBS -l nodes=$nnode:ppn=$ncore

cd \$PBS_O_WORKDIR/
module purge
module load g16-a03
module load anaconda3.7

start=\$(date "+%s.%N")

EOF
}

function main(){
start_job
cat << EOF >> mol_calc.sh
for molecule_name in \$(cat < molecule_name.txt)
do  
EOF
python_create_s0_gjf # create s0.gjf
g16job s0 opt-freq # calculate s0-opt-freq
python_mol_sub opt-freq # generate s1-opt-freq.gjf
g16job s1 opt-freq # calculate s1-opt-freq
get_trans_dip_adiabatic_energy #save as transition dipole, adiabatic energy as .npy file
python_mol_sub nacme # generate s1-nacme.gjf
g16job s1 nacme # calculate s1-nacme
dushin_job # generate dushin_cart.dus
nacme_job # generate nacme.out
# cp_to_ICrate_workdir #copy to IC calculation work dir
cat << EOF >> mol_calc.sh
done
EOF
}

main
qsub mol_calc.sh
