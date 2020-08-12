DIR=log_$(date +'%Y%m%d_%H%M%S')_$(hostname)/
# DIR=logtest_$(date +'%Y%m%d_%H%M%S')_$(hostname)/
mkdir -p $DIR
echo "Directory $DIR is created."
cp run.sh ${DIR}run.sh
 
cnt=0
for seed in 1
do
	# for id in 2 3 4 5 6 106 107 108 109 115 201 202
	# for id in 108 109 201 202
	for id in 204 205 206 207 208 209 210 211 212 214 215 216 217 218 219 220 221 224 225 227 228 232 233 234 235 236 237 239 240 241 242 245 246 247 249 250 252 253 254 256 258 259 260 261 263 264 265 266 267 268 271 272 273 274 275 276 277 278 279 282 283 284 285 288 291 292 293 296 297 299 300
	do
		for num_obs in 10
		do
			for noise_level in 0.0001
			do
				INS=i${id}_s${seed}_n${num_obs}_z${noise_level}/
				mkdir -p $DIR$INS

				for time_limit in 600 3600 10800
				do
					TIME=t${time_limit}/
					mkdir -p $DIR$INS$TIME

					model="heur"
					for init_solve in 1 #2
					do
						for stepsize in 3 #5
						do
							for fixlevel in 1 #2
							do
								ALGO=${model}_s${stepsize}_l${fixlevel}_i${init_solve}
								nohup julia main.jl $DIR$INS$TIME$ALGO $id $seed $num_obs $noise_level $time_limit $model $init_solve $stepsize $fixlevel > \
									$DIR$INS$TIME$ALGO.console 2>&1 &
								cnt=$(( cnt + 1 ))
								continue
							done
						done
					done

					model="minlp"
					for max_depth in 4 #2 3 4 5 6 7
					do
						ALGO=${model}_D${max_depth}
						nohup julia main.jl $DIR$INS$TIME$ALGO $id $seed $num_obs $noise_level $time_limit $model $max_depth > $DIR$INS$TIME$ALGO.console 2>&1 &
						cnt=$(( cnt + 1 ))
					done
				done
			done
		done
	done	
done 
echo "$cnt instances are running."
