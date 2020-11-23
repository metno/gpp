import time
import numpy as np
import argparse
import gridpp

def main():
    parser = argparse.ArgumentParser(description='Runs gridpp benchmarks for processing performance')
    parser.add_argument('-j', type=int, help='Run multiple cores', dest='num_cores')

    # if len(sys.argv) == 1:
    #     parser.print_help()
    #     sys.exit(1)

    args = parser.parse_args()

    input = dict()
    grids = dict()
    for i in [100, 200, 1000, 2000, 10000]:
        input[i] = np.random.rand(i, i)
    for i in [100, 200, 1000, 2000]:
        grids[i] = gridpp.Grid(*np.meshgrid(np.linspace(0, 1, i), np.linspace(0, 1, i)))
    radius = 7
    quantile = 0.5
    thresholds = np.linspace(0, 1, 11)
    run = dict()
    run[(gridpp.Grid, 1000)] = {"expected": 0.74, "args":np.meshgrid(np.linspace(0, 1, 1000), np.linspace(0, 1, 1000))}
    run[(gridpp.neighbourhood, 10000)] = {"expected": 2.27, "args":(input[10000], radius, gridpp.Mean)}
    run[(gridpp.neighbourhood,"2000 max")] = {"expected": 2.23, "args":(input[2000], radius, gridpp.Max)}
    run[(gridpp.neighbourhood_quantile_fast, 2000)] = {"expected": 1.41, "args":(input[2000], quantile, radius, thresholds)}
    run[(gridpp.bilinear, 200)] = {"expected": 1.43, "args":(grids[200], grids[200], input[200])}
    run[(gridpp.nearest, 1000)] = {"expected": 1.42, "args":(grids[1000], grids[1000], input[1000])}

    if args.num_cores is not None:
        print("Function                             Expected     Time     Diff    Scaling")
    else:
        print("Function                             Expected     Time     Diff")
    num_cores = [1]
    if args.num_cores is not None:
        num_cores += [args.num_cores]
    for key in run.keys()       :
        timings = dict()
        if isinstance(key, tuple):
            name = key[0].__name__ + " " + str(key[1])
            func = key[0]
        else:
            name = key.__name__
            func = key
        for num_core in num_cores:
            gridpp.set_omp_threads(num_core)
            s_time = time.time()
            func(*run[key]["args"])
            e_time = time.time()
            curr_time = e_time - s_time
            timings[num_core] = curr_time
            # print("%s() Expected: %.2f s Time: %.2f s" % (func.__name__, run[key]["expected"], e_time - s_time))
        diff = (timings[1] - run[key]["expected"]) / run[key]["expected"] * 100
        string = "%-36s %8.2f %8.2f %8.2f %%" % (name, run[key]["expected"], timings[1], diff)
        if args.num_cores is not None:
            scaling = timings[1] / timings[args.num_cores] / args.num_cores
            expected = timings[1] / args.num_cores
            scaling = 1 - (timings[args.num_cores] - expected) / (timings[1] - expected)
            # scaling = (1 - timings[args.num_cores] / timings[1]) * (args.num_cores + 1)

            string += " %8.2f %%" % (scaling * 100)
        print(string)
    # gridpp.neighbourhood(input, radius, gridpp.Mean)


if __name__ == "__main__":
    main()
