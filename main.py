import numpy as np
import matplotlib.pyplot as plt

def lotkaVolterra(t, x, p):
    # x[0] number of prey
    # x[1] number of predators
    a = p[0] # Reproductionrate of prey
    b = p[1] # Eatingrate of predators
    c = p[2] # Deathrate of predators
    d = p[3] # Reproductionrate of predators

    return np.array([ x[0] * (a - b * x[1]),
                     -x[1] * (c - d * x[0])])

def rkdp45(t0, y0, h, f, T, nmax, **kwargs):
    """
    Dormand-Prince method with adaptive stepsize control
    t0: initial time
    y0: initial value
    h: first timestep
    f: function
    T: end of time interval
    nmax: maximal amount of iterations
    tol: tolerance for truncation error
    """
    y = np.empty((nmax + 1, len(y0)), dtype = np.float64)
    y[0, :] = y0
    t = np.empty(nmax + 1, dtype = np.float64)
    t[0] = t0
    tol = kwargs.get("tolerance", 1e-14)
    parameters = kwargs.get("parameters",None)
    i = 0
    while i < nmax:
        k1 = f(t[i], y[i, :], parameters)
        k2 = f(t[i] + h / 5, y[i, :] + h * (1 / 5 * k1), parameters)
        k3 = f(t[i] + 3 / 10 * h, y[i, :] + h * (3 / 40 * k1 + 9 / 40 * k2), parameters)
        k4 = f(t[i] * 4 / 5 * h, y[i, :] + h * (44 / 45 * k1 - 56 / 15 * k2 + 32 / 9 * k3), parameters)
        k5 = f(t[i] * 8 / 9 * h, y[i, :] + h * (19372 / 6561 * k1 - 25360 / 2187 * k2 + 64448 / 6561 * k3 - 212 / 729 * k4), parameters)
        k6 = f(t[i] + h, y[i, :] + h * (9017 / 3168 * k1 - 355 / 33 * k2 + 46732 / 5247 * k3 + 49 / 176 * k4 - 5103 / 18656 * k5), parameters)
        u5 = y[i, :] + h * (35 / 384 * k1 + 500 / 1113 * k3 + 125 / 192 * k4 - 2187 / 6784 * k5 + 11 / 84 * k6)
        k7 = f(t[i] + h, u5, parameters)
        u4 = y[i, :] + h * (5179 / 57600 * k1 + 7571 / 16695 * k3 + 393 / 640 * k4 - 92097 / 339200 * k5 + 187 / 2100 * k6 + 1 / 40 * k7)

        # Truncationerror
        u_err = np.linalg.norm(u5-u4,np.inf)

        # Parameters for adaptive stepsize
        a_max = 2.0 # new step can be a max of two times as big
        a_min = 1 / 2 # new step can be min a half as big
        b = 0.8 # new time step will be multiplied by 0.8 to account for errors when calculating stepsize

        if u_err == 0:
            h_new = a_max * h
        else:
            h_new = h * min(a_max, max(a_min, b * np.power(tol / u_err, 1 / 5)))
        # if truncation error bigger than toleranc, try again
        if u_err > tol:
            h = h_new
            continue
        # last stepsize cannot overshoot T
        h = min(T - t[i], h_new)

        t[i + 1] = t[i] + h
        y[i + 1, :] = u5
        # If T is reached
        if T - t[i] <= 1e-30:
            return t[0:i + 1], y[0:i + 1, :]
        i += 1
    print("Max iterations reached at t =",t[-1])
    return t, y

if __name__ == "__main__":
    # Initial values
    t0 = 0
    T = 50

    prey = 5
    predators = 10

    # Parameters
    rp1 = 1.1 # Reproductionrate of Prey
    er = 0.4 # Eatingrate of Predators
    dr = 0.4 # Deathrate of Predators
    rp2 = 0.1 # Reproductionrte of Predators

    # Some other stuff
    y0 = np.array([prey,predators]) # Initial value
    h0 = 1e-5 # First timestep
    n = int((T-t0)/h0) * 2 # Max iterations

    # Dormand Prince Method
    t,y = rkdp45(t0,y0,h0,lotkaVolterra,T,n, tolerance = 1e-14, parameters = [rp1,er,dr,rp2])

    # Plot
    fig, ax = plt.subplots()
    ax.plot(t,y[:,0], label="Prey",color="tab:orange")
    ax.plot(t,y[:,1], label="Predator",color="tab:blue")
    ax.set_title("Lotka-Volterra\nParameters: a = " + str(rp1) + " b = " + str(er) + " c = " + str(dr) + " d = " + str(rp2))
    ax.set_xlabel("Time")
    ax.set_ylabel("Population")
    ax.legend()
    ax.grid()

    plt.show()