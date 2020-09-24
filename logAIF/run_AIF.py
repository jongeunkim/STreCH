import os
import sys
import time
sys.path.append('/home/ac.jnkim/development/AI-Feynman')

from feynman import run_aifeynman

easy_instances = [6, 8, 11, 12, 15, 16, 19, 20, 25, 27, 28, 35, 36, 37, 39, 40, 41, 46, 49, 50, 53, 54, 58, 59, 60, 65, 67, 68, 72, 73, 74, 75, 76, 78, 83, 85, 88, 92, 96, 97]
hard_instances = [4, 5, 7, 9, 10, 14, 17, 18, 21, 24, 32, 33, 34, 42, 45, 47, 52, 56, 61, 63, 64, 66, 71, 77, 79, 82, 84, 91, 93, 99, 100]
explog_instances = [1, 2, 3, 43, 44, 48, 80, 86, 87, 94]


def main():
    DIR = sys.argv[1]
    # id = int(sys.argv[2])
    # numobs = int(sys.argv[3])
    ops = sys.argv[2]
    # seed = int(sys.argv[2])
    # noise = float(sys.argv[4])

    # id=10
    # numobs=10
    # ops = "19ops"

    # INS = f"i{id}_n{numobs}_z0"
    # INS = f"i{id}_s{seed}_n{numobs}_z{noise}_{ops}"
    os.chdir(DIR)

    t_begin = time.time()
    run_aifeynman("", f"data.obs", 30, ops + ".txt", polyfit_deg=3, NN_epochs=500)
    t_elapsed = time.time() - t_begin

    with open(f"results/solution_data.obs", "r") as fp:
        first_line = fp.readline()
        for last_line in fp:
            pass

    ## Save the result
    with open('../result.txt', 'a+') as fp:
        fp.write(f"{DIR}\t{t_elapsed}\t{last_line}")

if __name__ == "__main__":
    main()