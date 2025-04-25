# 🚀 THERMAT – Thermal Analysis Tool for Rocket Engine Thrust Chambers

**THERMAT** is a Python-based tool designed to analyze heat transport in rocket engine combustion chambers and to support the design of regenerative cooling systems. 

It uses **Cantera** for chemical equilibrium and flow calculations, and **CoolProp** for accurate fluid property data.

---

## 📁 Project Structure

```plaintext
thermat/
├── main.py                  # Entry point: runs thermal analysis
├── core/              
│   ├── therm_functions.py	 # Main functions for geometry, flow, and thermal analysi
├── utilities/		           # Additional scripts not required for the core functionality
├── README.md
├── requirements.txt
```

---

## 🛫 Getting Started

### 1. Clone the repository

```bash
git clone https://github.com/flxbrt/thermat.git
cd thermat
```

### 2. Install dependencies

```bash
pip install -r requirements.txt
```

### 3. Run the tool
```bash
python main.py
```

This tool will
- Generate chamber and cooling channel geometry
- Solve flow properties 1D (axially discretized) based on **isentropic** flow along the combustion chamber; as of now the flow is computed for frozen conditions only; equilibirum flow will be implemented at a leter point
- Perform a thermal analysis and compute heat fluxes and temperatures (either fixed wall or regenerative cooling)
- Output plots and performance metrics such as absorbed heat and coolant heat capacity


---

## 💡 Features
- Allows fixed wall temperature or full regenerative cooling analysis with variable cooling geometry and number of channels
- Automatically detects infeasible cooling geometries
- Suports different propellant combinations
- Supports **frozen** flow / chemistry
- The hot gas side heat transfer is validated against the subscale SSME engine

---

## 🚀 Planned Enhancements
- **Equilibrium** flow / chemistry
- Provide calibration data on a 3kN liquid ethanol cryogenic oxygen engine once test campaigns has been conducted

---

## Author
Developed by @flxbrt
Version: 2.0
