## Set the time limit
TIMELIMIT=600

## Create a log folder
DIR=../log_$(date +'%Y%m%d_%H%M%S')_$(hostname)_t${TIMELIMIT}/
mkdir -p $DIR
echo "Directory $DIR is created."

## Copy the current codes
cp -R ../code/ ${DIR}/code/

## Set the parameter for data
id=9
seed=1
num_obs=10
noise_level=0.0001

## Create the counter for the instances
cnt=0

## Create a folder for a single data
INS=i${id}_s${seed}_n${num_obs}_z${noise_level}/
mkdir -p $DIR$INS

## Run heuristics
# method="heur"
# for max_depth in 3 #4 5
# do
# 	for init_solve in 2 # 1 2 3 4
# 	do
# 		ALGO=${method}_D${max_depth}_i${init_solve}
# 		nohup julia main.jl $DIR$INS$ALGO $id $seed $num_obs $noise_level $TIMELIMIT $method $max_depth $init_solve > $DIR$INS$ALGO.console 2>&1 &
# 		cnt=$(( cnt + 1 ))
# 		continue
# 	done
# done

# ## Run minlps
# method="minlp"
# for max_depth in 2 #3 4
# do
# 	for formulation in Cozad Cozad-CR New New-NR
# 	do
# 		ALGO=${method}_D${max_depth}_F-${formulation}
# 		nohup julia main.jl $DIR$INS$ALGO $id $seed $num_obs $noise_level $TIMELIMIT $method $max_depth $formulation > $DIR$INS$ALGO.console 2>&1 &
# 		cnt=$(( cnt + 1 ))
# 	done
# done

## Run optcheck
method="optcheck"
ALGO=${method}
nohup julia main.jl $DIR$INS$ALGO $id $seed $num_obs $noise_level $TIMELIMIT $method > $DIR$INS$ALGO.console 2>&1 &
cnt=$(( cnt + 1 ))

## Display how many instances are running
echo "$cnt instances are running."
