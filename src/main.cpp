#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <vector>
#include <getopt.h>

std::vector<double> inputToPhase(const std::vector<double> &data, double rate,
                                 const std::string &dataType);
std::vector<double> frequencyToPhase(const std::vector<double> &freqData,
                                     double rate);

std::tuple<double, double, int> calcAdevPhase(const std::vector<double> &phase,
                                              double rate, int mj, int stride);
void removeSmallNs(const std::vector<double> &taus,
                   const std::vector<double> &devs,
                   const std::vector<double> &deverrs,
                   const std::vector<int> &ns, std::vector<double> &outTaus,
                   std::vector<double> &outDevs,
                   std::vector<double> &outDeverrs, std::vector<int> &outNs);

std::tuple<std::vector<double>, std::vector<int>, std::vector<double>>
tauGenerator(const std::vector<double> &data, double rate,
             std::optional<std::string> tausMode = std::nullopt, bool v = false,
             bool even = false, int maximumM = -1) {

  if (rate == 0) {
    throw std::runtime_error("Warning! rate == 0");
  }

  std::vector<double> taus;
  if (!tausMode.has_value() || tausMode.value() == "octave") {
    int maxn = static_cast<int>(std::floor(std::log2(data.size())));
    for (int i = 0; i <= maxn; ++i) {
      taus.push_back(std::pow(2.0, i) / rate);
    }
  } else if (tausMode.value() == "log10") {
    double maxn = std::log10(data.size());
    int steps = static_cast<int>(10 * maxn);
    for (int i = 0; i < steps; ++i) {
      taus.push_back(std::pow(10.0, i / 10.0) / rate);
    }
    if (v) {
      std::cout << "tauGenerator: maxn=" << maxn << std::endl;
      for (auto tau : taus) {
        std::cout << " " << tau;
      }
      std::cout << std::endl;
    }
  } else if (tausMode.value() == "decade") {
    int maxn = static_cast<int>(std::floor(std::log10(data.size())));
    for (int k = 0; k <= maxn; ++k) {
      taus.push_back(1.0 * std::pow(10.0, k) / rate);
      taus.push_back(2.0 * std::pow(10.0, k) / rate);
      taus.push_back(4.0 * std::pow(10.0, k) / rate);
    }
  }

  std::vector<int> m;
  for (double tau : taus) {
    int mj = static_cast<int>(std::round(tau * rate));
    if (mj > 0 && mj < (maximumM == -1 ? data.size() : maximumM)) {
      m.push_back(mj);
    }
  }

  if (even) {
    m.erase(std::remove_if(m.begin(), m.end(),
                           [](int val) { return val % 2 != 0; }),
            m.end());
  }

  std::vector<double> tausUsed;
  for (int mj : m) {
    tausUsed.push_back(mj / rate);
  }

  if (v) {
    std::cout << "tauGenerator: m = ";
    for (int mj : m) {
      std::cout << mj << " ";
    }
    std::cout << std::endl;
  }

  return {data, m, tausUsed};
}

std::tuple< std::vector<double>, std::vector<double>, std::vector<double>, std::vector<int> >
adev(const std::vector<double> &data, double rate, const std::string &dataType,
     const std::optional<std::string> &tausMode = std::nullopt) {

  auto phase = inputToPhase(data, rate, dataType);

  auto [taus, m, tausUsed] = tauGenerator(phase, rate, tausMode);

  std::vector<double> ad(tausUsed.size(), 0.0);
  std::vector<double> ade(tausUsed.size(), 0.0);
  std::vector<int> adn(tausUsed.size(), 0);

  for (size_t idx = 0; idx < m.size(); ++idx) {
    try {
      auto [dev, deverr, n] = calcAdevPhase(phase, rate, m[idx], m[idx]);
      ad[idx] = dev;
      ade[idx] = deverr;
      adn[idx] = n;
    } catch (const std::runtime_error &e) {
      // std::cerr << "Error in calcAdevPhase for m[" << idx << "]: " << e.what() << std::endl;
      ad[idx] = 0.0;
      ade[idx] = 0.0;
      adn[idx] = 0;
    }
  }

  std::vector<double> outTaus, outDevs, outDeverrs;
  std::vector<int> outNs;
  removeSmallNs(tausUsed, ad, ade, adn,
                outTaus, outDevs,
                outDeverrs, outNs);

  return {outTaus, outDevs, outDeverrs, outNs};
}

std::tuple< std::vector<double>, std::vector<double>, std::vector<double>, std::vector<int> >
oadev(const std::vector<double> &data, double rate, const std::string &dataType,
     const std::optional<std::string> &tausMode = std::nullopt) {

  auto phase = inputToPhase(data, rate, dataType);

  auto [taus, m, tausUsed] = tauGenerator(phase, rate, tausMode);

  std::vector<double> ad(tausUsed.size(), 0.0);
  std::vector<double> ade(tausUsed.size(), 0.0);
  std::vector<int> adn(tausUsed.size(), 0);

  for (size_t idx = 0; idx < m.size(); ++idx) {
    try {
      auto [dev, deverr, n] = calcAdevPhase(phase, rate, m[idx], 1); // overlapping Allan Deviaton => stride=1
      ad[idx] = dev;
      ade[idx] = deverr;
      adn[idx] = n;
    } catch (const std::runtime_error &e) {
      // std::cerr << "Error in calcAdevPhase for m[" << idx << "]: " << e.what() << std::endl;
      ad[idx] = 0.0;
      ade[idx] = 0.0;
      adn[idx] = 0;
    }
  }

  std::vector<double> outTaus, outDevs, outDeverrs;
  std::vector<int> outNs;
  removeSmallNs(tausUsed, ad, ade, adn,
                outTaus, outDevs,
                outDeverrs, outNs);

  return {outTaus, outDevs, outDeverrs, outNs};
}
std::vector<double> inputToPhase(const std::vector<double> &data, double rate,
                                 const std::string &dataType) {
  if (dataType == "phase") {
    return data;
  } else if (dataType == "freq") {
    return frequencyToPhase(data, rate);
  } else {
    throw std::invalid_argument("Unknown data type: " + dataType);
  }
}

std::vector<double> frequencyToPhase(const std::vector<double> &freqData,
                                     double rate) {
  double dt = 1.0 / rate;
  std::vector<double> freqDataCopy = freqData;
  double mean = std::accumulate(freqDataCopy.begin(), freqDataCopy.end(), 0.0) /
                freqDataCopy.size();
  std::transform(freqDataCopy.begin(), freqDataCopy.end(), freqDataCopy.begin(),
                 [mean](double val) { return val - mean; });

  std::vector<double> phaseData(freqDataCopy.size() + 1, 0.0);
  for (size_t i = 1; i < phaseData.size(); ++i) {
    phaseData[i] = phaseData[i - 1] + freqDataCopy[i - 1] * dt;
  }
  return phaseData;
}

std::tuple<double, double, int> calcAdevPhase(const std::vector<double> &phase,
                                              double rate, int mj, int stride) {
  std::vector<double> d2, d1, d0;

  for (size_t i = 2 * mj; i < phase.size(); i += stride)
    d2.push_back(phase[i]);

  for (size_t i = mj; i < phase.size(); i += stride)
    d1.push_back(phase[i]);

  for (size_t i = 0; i < phase.size(); i += stride)
    d0.push_back(phase[i]);

  int n = std::min({d0.size(), d1.size(), d2.size()});

  if (n == 0) {
    throw std::runtime_error("Data array length is too small.");
  }

  double s = 0.0;
  for (int i = 0; i < n; ++i) {
    double v = d2[i] - 2 * d1[i] + d0[i];
    s += v * v;
  }

  double dev = std::sqrt(s / (2.0 * n)) / mj * rate;
  double deverr = dev / std::sqrt(n);

  return {dev, deverr, n};
}

void removeSmallNs(const std::vector<double> &taus,
                   const std::vector<double> &devs,
                   const std::vector<double> &deverrs,
                   const std::vector<int> &ns, std::vector<double> &outTaus,
                   std::vector<double> &outDevs,
                   std::vector<double> &outDeverrs, std::vector<int> &outNs) {
  for (size_t i = 0; i < ns.size(); ++i) {
    if (ns[i] > 1) {
      outTaus.push_back(taus[i]);
      outDevs.push_back(devs[i]);
      outDeverrs.push_back(deverrs[i]);
      outNs.push_back(ns[i]);
    }
  }
}

int main(int argc, char* argv[]) {
  std::string inputFileName;
  double samplePeriod = 0.0;
  std::string dataType;

  int opt;
  while ((opt = getopt(argc, argv, "i:s:t:")) != -1) {
    switch (opt) {
      case 'i':
        inputFileName = optarg;
        break;
      case 's':
        samplePeriod = std::stod(optarg);
        break;
      case 't':
        dataType = optarg;
        break;
      default:
        std::cerr << "Usage: " << argv[0] << " -i <inputFile> -s <samplePeriod> -t <dataType>" << std::endl;
        return 1;
    }
  }

  if (inputFileName.empty() || samplePeriod <= 0.0 || (dataType != "phase" && dataType != "freq")) {
    std::cerr << "Usage: " << argv[0] << " -i <inputFile> -s <samplePeriod> -t <dataType>" << std::endl;
    return 1;
  }

  std::cout << "# AllanToolsCxx - Overlapped Allan Deviaton test" << std::endl;
  std::cout << "# ==============================================" << std::endl;

  std::ifstream inputFile(inputFileName);
  if (!inputFile.is_open()) {
    throw std::runtime_error("Unable to open file");
  }

  std::vector<double> data;
  double value;
  while (inputFile >> value) {
    data.push_back(value);
  }

  auto [taus, adeviation, _, __] = oadev(data, 1.0 / samplePeriod, dataType);

  for (size_t i = 0; i < taus.size(); ++i) {
    std::cout << taus[i] << " " << adeviation[i] << std::endl;
  }

  return 0;
}
