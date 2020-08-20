# sleep 12000

TIMELIMIT=3600
# SLEEPTIME=10

DIR=../log_$(date +'%Y%m%d_%H%M%S')_$(hostname)_t${TIMELIMIT}/
mkdir -p $DIR
echo "Directory $DIR is created."
# cp run.sh ${DIR}run.sh
cp -R ../code/ ${DIR}/code/

cnt=0
for seed in 1
do
	# for id in 237
	# for id in 1
	## My own and from rational approximation paper, 10 instances
	# for id in 2 3 4 5 6 106 107 108 109 115
	## 40 instances, AIF with two or three variables
	# for id in 206 208 211 212 215 216 219 220 225 227 
	# for id in 228 235 236 237 239 240 241 246 249 250 
	# for id in 253 254 258 259 260 265 267 268 272 273
	# for id in 274 275 276 278 283 285 288 292 296 297
	## 31 instances, AIF with at least four variables
	for id in 204 205 207 209 210 214 217 218 221 224 232 233 234 242 245 247 252 256 261 263 264 266 271 277 279 282 284 291 293 299 300
	# for id in 204 205 207 209 210 214 217 218 221 224 232 233 234 242 245 247 
	# for id in 252 256 261 263 264 266 271 277 279 282 284 291 293 299 300
	## 10 instances, AIF with exp and/or log
	# for id in 201 202 203 243 244 248 280 286 287 294
	do
		for num_obs in 10
		do
			for noise_level in 0.0001
			do
				INS=i${id}_s${seed}_n${num_obs}_z${noise_level}/
				mkdir -p $DIR$INS

				model="heur"
				for max_depth in 3 4 5
				do
					for init_solve in 2 # 1 2 3 4
					do
						ALGO=${model}_D${max_depth}_i${init_solve}
						nohup julia main.jl $DIR$INS$ALGO $id $seed $num_obs $noise_level $TIMELIMIT $model $max_depth $init_solve > $DIR$INS$ALGO.console 2>&1 &
						cnt=$(( cnt + 1 ))
						continue
					done
				done

				model="minlp"
				for max_depth in 2 3 4
				do
					ALGO=${model}_D${max_depth}
					nohup julia main.jl $DIR$INS$ALGO $id $seed $num_obs $noise_level $TIMELIMIT $model $max_depth > $DIR$INS$ALGO.console 2>&1 &
					cnt=$(( cnt + 1 ))
				done

				# model="optcheck"
				# ALGO=${model}
				# nohup julia main.jl $DIR$INS$ALGO $id $seed $num_obs $noise_level $TIMELIMIT $model > $DIR$INS$ALGO.console 2>&1 &
				# cnt=$(( cnt + 1 ))
			done
		done
	done	
done 
echo "$cnt instances are running."
