#!/bin/bash
#PBS -N mol_sub_default
#PBS -o mol_sub_default.out
#PBS -e mol_sub_default.err
#PBS -q mem192v2
#PBS -l nodes=1:ppn=32

cd $PBS_O_WORKDIR/
module purge
module load g16-a03
module load anaconda3.7

start=$(date "+%s.%N")

for molecule_name in $(cat < molecule_name.txt)
do  
    python generate_gjf_from_name.py $molecule_name b3lyp &

    for job in $(jobs -p)
    do
        wait $job
    done

    INPUTFILE=$molecule_name/$molecule_name-b3lyp-s0-opt-freq.gjf

    # Setup G09 runtime environment
    OUTPUTFILE=${INPUTFILE%.*}.log
    export GAUSS_SCRDIR=/tmp/$USER/$RANDOM
    mkdir -p $GAUSS_SCRDIR

    echo "======================================================="
    echo "Working directory is" $PBS_O_WORKDIR
    echo "Starting on $(hostname) at $(date)"


    # Run g16 job
    cd $PBS_O_WORKDIR/
    time g16 < $INPUTFILE > $OUTPUTFILE

    # Clean scratch file once finished
    rm -rf $GAUSS_SCRDIR

    for job in $(jobs -p)
    do
        wait $job
    done

    # Check whether process error
    check=$(tail -n 1 $molecule_name/$molecule_name-b3lyp-s0-opt-freq.log | grep 'Normal termination')
    if [ -z "$check" ]; then 
      echo "s0-opt-freq process is not normally terminated."
      continue
    fi

    formchk $molecule_name-b3lyp-s0-opt-freq.chk
    mv $molecule_name-b3lyp-s0-opt-freq* $molecule_name/

    python runscript_mol_sub.py $molecule_name opt-freq b3lyp&

    for job in $(jobs -p)
    do
        wait $job
    done
    INPUTFILE=$molecule_name/$molecule_name-b3lyp-s1-opt-freq.gjf

    # Setup G09 runtime environment
    OUTPUTFILE=${INPUTFILE%.*}.log
    export GAUSS_SCRDIR=/tmp/$USER/$RANDOM
    mkdir -p $GAUSS_SCRDIR

    echo "======================================================="
    echo "Working directory is" $PBS_O_WORKDIR
    echo "Starting on $(hostname) at $(date)"


    # Run g16 job
    cd $PBS_O_WORKDIR/
    time g16 < $INPUTFILE > $OUTPUTFILE

    # Clean scratch file once finished
    rm -rf $GAUSS_SCRDIR

    for job in $(jobs -p)
    do
        wait $job
    done

    # Check whether process error
    check=$(tail -n 1 $molecule_name/$molecule_name-b3lyp-s1-opt-freq.log | grep 'Normal termination')
    if [ -z "$check" ]; then 
      echo "s1-opt-freq process is not normally terminated."
      continue
    fi

    formchk $molecule_name-b3lyp-s1-opt-freq.chk
    mv $molecule_name-b3lyp-s1-opt-freq* $molecule_name/

    python runscript_get_dipole_and_energy.py $molecule_name b3lyp
    for job in $(jobs -p)
    do
        wait $job
    done

    sleep 1
    python runscript_mol_sub.py $molecule_name nacme b3lyp&

    for job in $(jobs -p)
    do
        wait $job
    done
    INPUTFILE=$molecule_name/$molecule_name-b3lyp-s1-nacme.gjf

    # Setup G09 runtime environment
    OUTPUTFILE=${INPUTFILE%.*}.log
    export GAUSS_SCRDIR=/tmp/$USER/$RANDOM
    mkdir -p $GAUSS_SCRDIR

    echo "======================================================="
    echo "Working directory is" $PBS_O_WORKDIR
    echo "Starting on $(hostname) at $(date)"


    # Run g16 job
    cd $PBS_O_WORKDIR/
    time g16 < $INPUTFILE > $OUTPUTFILE

    # Clean scratch file once finished
    rm -rf $GAUSS_SCRDIR

    for job in $(jobs -p)
    do
        wait $job
    done

    # Check whether process error
    check=$(tail -n 1 $molecule_name/$molecule_name-b3lyp-s1-nacme.log | grep 'Normal termination')
    if [ -z "$check" ]; then 
      echo "s1-nacme process is not normally terminated."
      continue
    fi

    formchk $molecule_name-b3lyp-s1-nacme.chk
    mv $molecule_name-b3lyp-s1-nacme* $molecule_name/

    cp $PBS_O_WORKDIR/$molecule_name/$molecule_name-b3lyp-*-opt-freq.fchk /home/schsu/MOMAP_bak/momap/dushin
    cp $PBS_O_WORKDIR/$molecule_name/$molecule_name-b3lyp-*-opt-freq.log /home/schsu/MOMAP_bak/momap/dushin

    cd /home/schsu/MOMAP_bak/momap/dushin
    sed -i '1s/g16/g09/g' /home/schsu/MOMAP_bak/momap/dushin/$molecule_name-b3lyp-s0-opt-freq.log
    sed -i '1s/g16/g09/g' /home/schsu/MOMAP_bak/momap/dushin/$molecule_name-b3lyp-s1-opt-freq.log

    change8="2 1 'ground state' '$molecule_name-b3lyp-s1-opt-freq.log'"
    change9="1 1 'excited state' '$molecule_name-b3lyp-s0-opt-freq.log'"
    sed -i "8c $change8" /home/schsu/MOMAP_bak/momap/dushin/dushin.inp
    sed -i "9c $change9" /home/schsu/MOMAP_bak/momap/dushin/dushin.inp

    source /home/schsu/MOMAP_bak/momap/dushin/run 2> $PBS_O_WORKDIR/$molecule_name/dushin_res.txt

    cp /home/schsu/MOMAP_bak/momap/dushin/dushin_cart.dus $PBS_O_WORKDIR/$molecule_name/

    rm /home/schsu/MOMAP_bak/momap/dushin/$molecule_name*
    rm dushin_cart.dus
    cp $PBS_O_WORKDIR/$molecule_name/$molecule_name-b3lyp-s1-opt-freq.log /home/schsu/MOMAP_bak/momap/nacme
    cp $PBS_O_WORKDIR/$molecule_name/$molecule_name-b3lyp-s1-nacme.log /home/schsu/MOMAP_bak/momap/nacme

    cd /home/schsu/MOMAP_bak/momap/nacme
    sed -i "5s/get-nacme.*/get-nacme $molecule_name-b3lyp-s1-nacme.log $molecule_name-b3lyp-s1-opt-freq.log/g" /home/schsu/MOMAP_bak/momap/nacme/run

    source /home/schsu/MOMAP_bak/momap/nacme/run 2> $PBS_O_WORKDIR/$molecule_name/nacme_res.txt
    cp /home/schsu/MOMAP_bak/momap/nacme/nacme.out $PBS_O_WORKDIR/$molecule_name/
    cp /home/schsu/MOMAP_bak/momap/nacme/nacme.log $PBS_O_WORKDIR/$molecule_name/

    rm /home/schsu/MOMAP_bak/momap/nacme/$molecule_name*
    rm nacme.out
    rm nacme.log
    
    for job in $(jobs -p)
    do
        wait $job
    done

    end=$(date "+%s.%N")
    runtime=$(echo "$end - $start" | bc -lq)

    echo "Job Ended at $(date)"
    echo "Estimated Runtime: $runtime seconds."
    echo '======================================================='
    #mpdallexit
done
