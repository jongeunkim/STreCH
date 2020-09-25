maxCPU=7
CPUid=1
# numlines=`wc -l < result.txt`
echo "Start the first batch"

seed=1
noise=0
# for id in 6 8 11 12 15 16 19 20 25 27 28 35 36 37 39 40 41 46 49 50
# for id in 53 54 58 59 60 65 67 68 72 73 74 75 76 78 83 85 88 92 96 97 
# for id in 4 5 
# for id in 7 9 10 14 17 18 21 24 32 33 34 42 45 47  
for id in 52 56 61 63 64 66 71 77 79 82 84 91 93 99 100
do
    for numobs in 10 #100 1000
    do
        for ops in 19ops #14ops
        do
            ## Create a folder
            DIR=i${id}_s${seed}_n${numobs}_z${noise}_${ops}/
            mkdir -p ${DIR}

            ## Generate a dataset
            julia obsgen.jl $DIR $id $seed $numobs $noise

            ## Run AIF
            nohup python run_AIF.py $DIR $ops > ${DIR}console.console 2>&1 &
            # echo "id=$id at cpu $CPUid"
            # nohup taskset --cpu-list ${CPUid} python test_AIF.py ${id} ${numobs} ${ops} > ${DIR}console.console 2>&1 &
            # CPUid=$((CPUid+1))

            # if [ $CPUid -gt $maxCPU ]
            # then
            #     numlines2=`wc -l < result.txt`
            #     while [ $numlines2 -lt $((numlines + maxCPU)) ]
            #     do
            #         echo "Wait for one minute" $numlines2 $((numlines + maxCPU))
            #         sleep 60
            #     done

            #     echo "Start the next batch"
            #     CPUid=1
            #     numlines=`wc -l < result.txt`
            # fi
        done
    done

done
