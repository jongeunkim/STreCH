DIR=$(date +'%Y%m%d_%H%M%S')/
mkdir $DIR
echo "Directory $DIR is created."

cnt=0
for seed in 1 2 3
do
	for id in 1 2 3 4 5 6
	do
		for num_obs in 10
		do
			for noise_level in 0.0001
			do
				for time_limit in 300 3600 10800
				do
					model="heur"
					max_depth=10
					for stepsize in 1 2 3
					do
						NAME=i${id}_s${seed}_n${num_obs}_z${noise_level}_D${max_depth}_t${time_limit}_${model}${stepsize}
						nohup julia main.jl $DIR$NAME $id $seed $num_obs $noise_level $max_depth $time_limit $model $stepsize > $DIR$NAME.console 2>&1 &
						cnt=$(( cnt + 1 ))
					done

					# model="minlp"
					# for max_depth in 0 1 2 3 4 5 
					# do
					# 	NAME=i${id}_s${seed}_n${num_obs}_z${noise_level}_D${max_depth}_t${time_limit}_${model}
					# 	nohup julia main.jl $DIR$NAME $id $seed $num_obs $noise_level $max_depth $time_limit $model 0 > $DIR$NAME.console 2>&1 &
					# 	cnt=$(( cnt + 1 ))
                    # done
				done
			done
		done
	done	
done 
echo "$cnt instances are running."
