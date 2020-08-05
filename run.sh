DIR=log_$(date +'%Y%m%d_%H%M%S')/
mkdir -p $DIR
echo "Directory $DIR is created."
cp run.sh ${DIR}run.sh
 
cnt=0
for seed in 1
do
	#for id in 1 2 3 4 5 6
	for id in 106 107 108 109 115
	do
		for num_obs in 10
		do
			for noise_level in 0.0001
			do
				INS=i${id}_s${seed}_n${num_obs}_z${noise_level}/
				mkdir -p $DIR$INS

				for time_limit in 300 3600 7200 10800
				do
					TIME=t${time_limit}/
					mkdir -p $DIR$INS$TIME

					model="heur"
					max_depth=10
					for stepsize in 3 4 5
					do
						ALGO=${model}_s${stepsize}_D${max_depth}
						nohup julia main.jl $DIR$INS$TIME$ALGO $id $seed $num_obs $noise_level $max_depth $time_limit $model $stepsize > $DIR$INS$TIME$ALGO.console 2>&1 &
						cnt=$(( cnt + 1 ))
					done

					model="minlp"
					for max_depth in 0 1 2 3 4 5 6 7
					do
						ALGO=${model}_D${max_depth}
						nohup julia main.jl $DIR$INS$TIME$ALGO $id $seed $num_obs $noise_level $max_depth $time_limit $model 0 > $DIR$INS$TIME$ALGO.console 2>&1 &
						cnt=$(( cnt + 1 ))
                    done
				done
			done
		done
	done	
done 
echo "$cnt instances are running."
