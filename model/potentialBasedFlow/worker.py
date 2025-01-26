
import random
import time


from amplpy import AMPL


# Define constants
# Economic parameters
cost = 2
retail = 15
recover = -3
# Parameters to generate data in each process
samples = 10000
sigma = 100
mu = 400

def worker(data):
    """
    Example worker
    Input: a list with parameters
        alpha - parameter to be passed to AMPL
        beta - parameter to be passed to AMPL
        run - number of the run for the alpha and beta combination
        seed - seed to generate data for the given process
    Output: a list with parameters
        alpha - parameter used in the run
        beta - parameter used in the run
        run - number of the run for the alpha and beta combination
        obj - objective value of the solved model
        worker_time - wall time used by the worker
    This function as the following steps: 
        generate data acording to the received parameters
        instantiate AMPL and load an existing model file
        solve the model with the defined solver
        output results
    """
    start_time = time.time()

    # Get information from the input list
    alpha = data[0]
    beta = data[1]
    run = data[2]
    seed = data[3]

    # Initialyze and seed random number generator
    rng = random.Random()
    rng.seed(seed)

    # Generate data for this execution
    demand = [max(rng.normalvariate(mu,sigma), 0) for i in range(samples)]
    maxrev = max(demand) * (retail - cost)
    minrev = max(demand) * (recover - cost) + min(demand) * retail

    # Create AMPL instance and load the model
    ampl = AMPL()
    ampl.read("model.mod")

    # Load the data
    ampl.param["samples"] = samples
    ampl.param["cost"] = cost
    ampl.param["recover"] = recover
    ampl.param["retail"] = retail
    ampl.param["minrev"] = minrev
    ampl.param["maxrev"] = maxrev
    ampl.param["demand"] = demand
    ampl.param["alpha"] = alpha
    ampl.param["beta"] = beta

    # Set solver and options
    ampl.option["solver"] = "highs"
    ampl.option["highs_options"] = "tech:threads=1"
    # Solve without generating any output
    # use this when your model is ready for deployment
    # otherwise use the regular solve (commented below) to track potential issues
    ampl.get_output("solve;")
    #ampl.solve()

    # Get results
    obj = ampl.obj["prof"].value()

    worker_time = time.time() - start_time
    return [alpha, beta, run, obj, worker_time]
