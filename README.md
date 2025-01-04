# AllanToolsCxx

## Overview
AllanToolsCxx is a C++ tool designed for calculating and analyzing Allan Deviation, which is widely used in frequency stability analysis. The project includes functionality to compute Allan Deviation and generate graphical visualizations, along with Python scripts for additional data processing and plotting.

## Features
- Compute Allan Deviation and Overlapping Allan Deviation.
- Command-line interface for easy integration into workflows.
- Python scripts for generating and visualizing Allan Deviation plots.
- Simple and modular design for extensibility.

## Project Structure
```plaintext
.
├── CMakeLists.txt          # Build configuration for the project
├── run.sh                 # Script to execute the Allan Deviation calculation
├── run_generate.sh        # Script to generate Allan Deviation results and save to a file
├── allan_plot.py          # Python script for plotting Allan Deviation results
├── show_allan_plot_test.sh # Script to display a test plot using the generated Allan Deviation data
├── allantools_test.py     # Python implementation of Allan Deviation calculation for comparison
├── src
│   └── main.cpp           # Main C++ implementation
```

## Getting Started

### Prerequisites
- C++17 compatible compiler (e.g., GCC, Clang).
- CMake version 3.10 or higher.
- Python 3 with `matplotlib` installed for plotting.

### Build Instructions
1. Clone the repository:
   ```bash
   git clone https://github.com/peterkaczorowski/AllanToolsCxx.git
   cd AllanToolsCxx
   ```

2. Create a build directory and configure the project:
   ```bash
   mkdir build && cd build
   cmake ..
   ```

3. Compile the project:
   ```bash
   make
   ```

4. Run the executable:
   ```bash
   ./AllanDeviationLab -i <input_file> -s <sample_period> -t <data_type>
   ```

### Example Usage
Run Allan Deviation calculation using a test dataset:
```bash
./run.sh
```

Generate Allan Deviation results and save them to a file:
```bash
./run_generate.sh
```

Visualize the results with a Python script:
```bash
./show_allan_plot_test.sh
```

### Python Plotting
Use the `allan_plot.py` script to plot Allan Deviation data:
```bash
python3 allan_plot.py -i <data_file>
```

## Author
Piotr Kaczorowski

## License
This project is licensed under the [MIT License](LICENSE).

## Acknowledgments
This project was inspired by the need for accurate frequency stability analysis tools and leverages Allan Deviation as a core metric for analysis.
The code is based on the AllanTools for Python project by Anders Wallin, which is available on GitHub at: https://github.com/aewallin/allantools.

