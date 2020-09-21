TIMELIMIT=10800			# time limit for a method (MINLP or Heuristic)
NUM_ID_PER_BATCH=18		# the number of instances in a batch. If you run multiple methods, then it runs "NUM_ID_PER_BATCH x NUM_METHODS" per batch. 

## Initializeation
cnt=0					# the number of runs, +1 when it runs a instance (=data+method)
cntid=0					# the number of ids, +1 when it runs all methods for a data 
sleeptime=0				# sleep time, we increase this per batch so that different batches do not run at the same time.

## Create a log folder
DIR=../log_$(date +'%Y%m%d_%H%M%S')_$(hostname)_t${TIMELIMIT}/
mkdir -p $DIR
echo "Directory $DIR is created."

## Copy the current codes
cp -R ../code/ ${DIR}/code/

for seed in 1
do
	# for id in 37
	for id in 9

	## easy40 + hard31
	# for id in 6 8 11 12 15 16 19 20 25 27 28 35 36 37 39 40 41 46 49 50 53 54 58 59 60 65 67 68 72 73 74 75 76 78 83 85 88 92 96 97 4 5 7 9 10 14 17 18 21 24 32 33 34 42 45 47 52 56 61 63 64 66 71 77 79 82 84 91 93 99 100

	## 40 instances, AIF with two or three variables
	# for id in 6 8 11 12 15 16 19 20 25 27 28 35 36 37 39 40 41 46 49 50 53 54 58 59 60 65 67 68 72 73 74 75 76 78 83 85 88 92 96 97
	# for id in 6 8 11 12 15 16 19 20 25 27  28 35 36 37 39 40 41 46 49 50 
	# for id in 53 54 58 59 60 65 67 68 72 73  74 75 76 78 83 85 88 92 96 97
	
	## 31 instances, AIF with at least four variables
	# for id in 4 5 7 9 10 14 17 18 21 24 32 33 34 42 45 47 52 56 61 63 64 66 71 77 79 82 84 91 93 99 100
	# for id in 4 5 7 9 10 14 17 18 21 24 32 33 34 42 45 47 
	# for id in 52 56 61 63 64 66 71 77 79 82 84 91 93 99 100

	## 10 instances, AIF with exp and/or log
	# for id in 1 2 3 43 44 48 80 86 87 94
	do
		for num_obs in 10
		do
			for noise_level in 0.0001
			do
				INS=i${id}_s${seed}_n${num_obs}_z${noise_level}/
				mkdir -p $DIR$INS

				# method="heur"
				# for max_depth in 3 4 5
				# do
				# 	for init_solve in 2 # 1 2 3 4
				# 	do
				# 		ALGO=${method}_D${max_depth}_i${init_solve}
				# 		nohup julia main.jl $DIR$INS$ALGO $id $seed $num_obs $noise_level $TIMELIMIT $method $max_depth $init_solve > $DIR$INS$ALGO.console 2>&1 &
				# 		cnt=$(( cnt + 1 ))
				# 		continue
				# 	done
				# done

				method="minlp"
				for max_depth in 2 #3 4
				do
					for formulation in Cozad Cozad-CR New New-NR
					do
						ALGO=${method}_D${max_depth}_F-${formulation}
						nohupout="$DIR$INS$ALGO.console"
						nohup sleep $sleeptime > $nohupout 2>&1 && 
							julia main.jl $DIR$INS$ALGO $id $seed $num_obs $noise_level $TIMELIMIT $method $max_depth $formulation > $nohupout 2>&1 &
						cnt=$(( $cnt + 1 ))
					done
				done

				# method="optcheck"
				# ALGO=${method}
				# nohup julia main.jl $DIR$INS$ALGO $id $seed $num_obs $noise_level $TIMELIMIT $method > $DIR$INS$ALGO.console 2>&1 &
				# cnt=$(( cnt + 1 ))
			done
		done

		cntid=$(( $cntid + 1))
		if (( $cntid % $NUM_ID_PER_BATCH  == 0 ))
		then
			echo "$cnt instances are running."
			sleeptime=$(( $sleeptime + $TIMELIMIT + 120 ))
		fi
	done
done 

echo "$cnt instances are running."
sleeptime=$(( $sleeptime + $TIMELIMIT ))
# echo "All runs will end in $(( $sleeptime / 3600 )) hours."
echo "All runs will end in $(( $sleeptime / (24 * 3600) )) days $(( $sleeptime / 3600 - $sleeptime / (24 * 3600) * 24 )) hours $(( $sleeptime / 60 - $sleeptime / 3600 * 60 )) minutes."