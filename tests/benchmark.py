import time
import numpy as np
import argparse
import gridpp

def main():
    parser = argparse.ArgumentParser(description='Runs gridpp benchmarks for processing performance')
    parser.add_argument('-j', type=int, help='Run multiple cores', dest='num_cores')
    parser.add_argument('-s', type=int, default=1, help='Enlarge the inputs by this scaling factor to run a bigger test', dest='scaling')
    parser.add_argument('-n', type=int, default=1, help='Number of iterations to average over', dest='iterations')
    parser.add_argument('-t', help='Run only this function', dest="function")

    # if len(sys.argv) == 1:
    #     parser.print_help()
    #     sys.exit(1)

    args = parser.parse_args()

    input = dict()
    grids = dict()
    points = dict()
    np.random.seed(1000)

    for i in [10, 50, 100, 200, 1000, 2000, 10000]:
        input[i] = np.random.rand(i * args.scaling, i)*10
    for i in [10, 50, 100, 200, 1000]:
        grids[i] = gridpp.Grid(*np.meshgrid(np.linspace(0, 1, i), np.linspace(0, 1, i * args.scaling)))
    for i in [1000]:
        # points[i] = gridpp.Points(np.linspace(0, 1, i), np.zeros(i))
        points[i] = gridpp.Points(np.random.rand(i) * 10, np.random.rand(i) * 10)
    structure = gridpp.BarnesStructure(10000)
    radius = 7
    quantile = 0.5
    thresholds = np.linspace(0, 1, 11)
    run = dict()
    run[(gridpp.Grid, "1000²")] = {"expected": 0.74, "args":np.meshgrid(np.linspace(0, 1, 1000 * args.scaling), np.linspace(0, 1, 1000))}
    run[(gridpp.neighbourhood, "10000²")] = {"expected": 2.05, "args":(np.zeros([10000, 10000]), radius, gridpp.Mean)}
    run[(gridpp.neighbourhood,"2000² max")] = {"expected": 2.23, "args":(input[2000], radius, gridpp.Max)}
    run[(gridpp.neighbourhood_quantile_fast, "2000²")] = {"expected": 1.23, "args":(input[2000], quantile, radius, thresholds)}
    run[(gridpp.bilinear, "1000²")] = {"expected": 1.68, "args":(grids[1000], grids[1000], input[1000])}
    run[(gridpp.nearest, "1000²")] = {"expected": 1.52, "args":(grids[1000], grids[1000], input[1000])}
    run[(gridpp.optimal_interpolation, "1000² 1000")] = {"expected": 1.57, "args":(grids[1000],
        input[1000], points[1000], np.zeros(1000), np.ones(1000), np.ones(1000), structure, 20)}
    run[(gridpp.dewpoint, "1e7")] = {"expected": 0.53, "args":(np.zeros(10000000) + 273.15, np.zeros(10000000))}

    print("Gridpp version %s" % gridpp.version())
    if args.num_cores is not None:
        print("Function                             Expected     Time     Diff    Scaling")
    else:
        print("Function                             Expected     Time     Diff")
    num_cores = [1]
    if args.num_cores is not None:
        num_cores += [args.num_cores]
    for key in run.keys()       :
        timings = dict()
        for num_core in num_cores:
            timings[num_core] = 0

        if isinstance(key, tuple):
            name = key[0].__name__ + " " + str(key[1])
            func = key[0]
        else:
            name = key.__name__
            func = key
        if args.function is not None:
            if func.__name__ != args.function:
                continue
        for num_core in num_cores:
            gridpp.set_omp_threads(num_core)
            for it in range(args.iterations):
                s_time = time.time()
                func(*run[key]["args"])
                e_time = time.time()
                curr_time = e_time - s_time
                timings[num_core] += curr_time
                # print("%s() Expected: %.2f s Time: %.2f s" % (func.__name__, run[key]["expected"], e_time - s_time))
        for num_core in num_cores:
            timings[num_core] /= args.iterations

        diff = (timings[1] - run[key]["expected"] * args.scaling) / (run[key]["expected"]  * args.scaling) * 100
        string = "%-36s %8.2f %8.2f %8.2f %%" % (name, run[key]["expected"] * args.scaling, timings[1], diff)
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
