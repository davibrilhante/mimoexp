DETECTORS=('zf' 'mmse' 'zfsic')
LENGTH=(256 512 1024)
ANTENNAS=(2 4 8)
ORDER=(4 16)
SNR=`seq -10 30`

for d in ${DETECTORS[@]}
do
    for l in ${LENGTH[@]}
    do
        for a in ${ANTENNAS[@]}
        do
            for o in ${ORDER[@]}
            do
                for s in ${SNR[@]}
                do
                    echo -e "python3 mqam_example.py --$d -l $l --tx $a --rx $a -m $o --snr $s >> $d-$l-$a-$o" >> submit
                done
            done
        done
    done
done
