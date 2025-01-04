import numpy as np


def adev(data, rate, data_type, taus=None):
    phase = input_to_phase(data, rate, data_type)
    (phase, m, taus_used) = tau_generator(phase, rate, taus)

    ad = np.zeros_like(taus_used)
    ade = np.zeros_like(taus_used)
    adn = np.zeros_like(taus_used)


    for idx, mj in enumerate(m):
        (ad[idx], ade[idx], adn[idx]) = calc_adev_phase(phase, rate, mj, mj)
        # print(f"{ad[idx]} | {ade[idx]} | {adn[idx]}")

    return remove_small_ns(taus_used, ad, ade, adn)

def calc_adev_phase(phase, rate, mj, stride):
    mj = int(mj)
    stride = int(stride)
    d2 = phase[2 * mj::stride]
    d1 = phase[1 * mj::stride]
    d0 = phase[::stride]

    n = min(len(d0), len(d1), len(d2))

    if n == 0:
        RuntimeWarning("Data array length is too small: %i" % len(phase))
        print(f"Data array length is too small: %i" % len(phase))
        n = 1

    v_arr = d2[:n] - 2 * d1[:n] + d0[:n]
    s = np.sum(v_arr * v_arr)

    dev = np.sqrt(s / (2.0*n)) / mj*rate
    deverr = dev / np.sqrt(n)

    # std::cout << "s:" << s << " n:" << n << " mj:" << mj << " rate:" << rate << std::endl;
    # print(f"s: {s} n:{n} mj:{mj} rate:{rate} => dev:{dev}")

    return dev, deverr, n

def input_to_phase(data, rate, data_type):
    if data_type == "phase":
        return data
    elif data_type == "freq":
        return frequency2phase(data, rate)
    else:
        raise Exception("unknown data_type: " + data_type)

def frequency2phase(freqdata, rate):
    dt = 1.0 / float(rate)
    freqdata = trim_data(freqdata)
    freqdata = freqdata - np.nanmean(freqdata)
    phasedata = np.cumsum(freqdata) * dt
    phasedata = np.insert(phasedata, 0, 0)
    return phasedata

def tau_generator(data, rate, taus=None, v=False, even=False, maximum_m=-1):
    if rate == 0:
        raise RuntimeError("Warning! rate==0")

    if taus is None:  # empty or no tau-list supplied
        taus = "octave"  # default to octave
    elif isinstance(taus, list) and taus == []:  # empty list
        taus = "octave"

    if isinstance(taus, np.ndarray) or isinstance(taus, list) and len(taus):
        pass
    elif taus == "all":  # was 'is'
        taus = (1.0/rate)*np.linspace(1.0, len(data), len(data))
    elif taus == "octave":
        maxn = np.floor(np.log2(len(data)))
        taus = (1.0/rate)*np.logspace(0, int(maxn), int(maxn+1), base=2.0)
    elif taus == "log10":
        maxn = np.log10(len(data))
        taus = (1.0/rate)*np.logspace(0, maxn, int(10*maxn), base=10.0)
        if v:
            print("tau_generator: maxn %.1f" % maxn)
            print("tau_generator: taus=" % taus)
    elif taus == "decade":  # 1, 2, 4, 10, 20, 40, spacing similar to Stable32
        maxn = np.floor(np.log10(len(data)))
        taus = []
        for k in range(int(maxn+1)):
            taus.append(1.0*(1.0/rate)*pow(10.0, k))
            taus.append(2.0*(1.0/rate)*pow(10.0, k))
            taus.append(4.0*(1.0/rate)*pow(10.0, k))

    data, taus = np.array(data), np.array(taus)
    rate = float(rate)
    m = []  # integer averaging factor. tau = m*tau0

    if maximum_m == -1:  # if no limit given
        maximum_m = len(data)

    m = np.round(taus * rate)
    taus_valid1 = m < len(data)
    taus_valid2 = m > 0
    taus_valid3 = m <= maximum_m
    taus_valid = taus_valid1 & taus_valid2 & taus_valid3
    m = m[taus_valid]
    m = m[m != 0]       # m is tau in units of datapoints
    m = np.unique(m)    # remove duplicates and sort

    if v:
        print("tau_generator: ", m)

    if len(m) == 0:
        print("Warning: sanity-check on tau failed!")
        print("   len(data)=", len(data), " rate=", rate, "taus= ", taus)

    taus2 = m / float(rate)

    if even:  # used by Theo1
        m_even_mask = ((m % 2) == 0)
        m = m[m_even_mask]
        taus2 = taus2[m_even_mask]

    return data, m, taus2

def remove_small_ns(taus, devs, deverrs, ns):
    ns_big_enough = ns > 1

    o_taus = taus[ns_big_enough]
    o_devs = devs[ns_big_enough]
    o_ns = ns[ns_big_enough]
    if isinstance(deverrs, list):
        assert len(deverrs) < 3
        o_deverrs = [deverrs[0][ns_big_enough], deverrs[1][ns_big_enough]]
    else:
        o_deverrs = deverrs[ns_big_enough]
    if len(o_devs) == 0:
        print("remove_small_ns() nothing remains!?")
        raise UserWarning

    return o_taus, o_devs, o_deverrs, o_ns



def main():
    print("allan deviation lab")
    print("===================")

    # set global options
    np.set_printoptions(threshold=np.inf)

    # load data
    data = np.loadtxt("HP866A-16Mhz.tic")

    # set data opts
    sample_period = 0.004735426008968611
    data_type = "phase"

    # calculate Allan Deviation
    taus, adeviation, _, _ = adev(data, rate=(1/sample_period), data_type=data_type)

    # print results
    for tau, deviation in zip(taus, adeviation):
        print(f"{tau:>10.6f} {deviation:>15.6e}" )


if __name__ == "__main__":
    main()
