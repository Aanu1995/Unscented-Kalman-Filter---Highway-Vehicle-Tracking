# Unscented Kalman Filter - Highway Vehicle Tracking

A sensor fusion implementation using Unscented Kalman Filter to track multiple vehicles on a highway with LiDAR and Radar measurements.

![UKF Highway Vehicle Tracking](media/ukf_highway_tracked.gif)

## Project Overview

This project implements an **Unscented Kalman Filter (UKF)** to track multiple vehicles in a simulated highway environment. The filter fuses noisy LiDAR and Radar sensor data to estimate vehicle positions and velocities with high accuracy.

### Key Features

- Multi-object tracking of 3 highway vehicles
- Sensor fusion combining LiDAR (position) and Radar (position + velocity)
- Real-time 3D visualization using PCL
- CTRV (Constant Turn Rate and Velocity) motion model
- Achieved RMSE: [0.061, 0.152, 0.392, 0.535] for [px, py, vx, vy]

## Implementation Methodology

### 1. Motion Model

- **State Vector**: 5D CTRV model [px, py, velocity, yaw_angle, yaw_rate]
- **Process Model**: Handles non-linear vehicle motion including turning
- **Process Noise**: Tuned parameters for longitudinal (2.4 m/s²) and yaw acceleration (2.9 rad/s²)

### 2. Unscented Kalman Filter Algorithm

**Prediction Step:**

- Generate 15 sigma points using augmented state (7D with process noise)
- Propagate sigma points through non-linear motion model
- Calculate predicted state mean and covariance

**Update Step:**

- **LiDAR**: Linear measurement model for [px, py] positions
- **Radar**: Non-linear measurement model for [range, bearing, range_rate]
- Separate update functions optimized for each sensor type

### 3. Sensor Fusion Strategy

- **LiDAR**: Provides precise position measurements with low noise
- **Radar**: Adds velocity information and works in all weather conditions
- **Combined**: Complementary sensors improve overall tracking accuracy

### 4. Parameter Optimization

Systematic tuning of process noise parameters through:

- Testing 17 different parameter combinations
- Minimizing RMSE across all tracked vehicles
- Validating filter consistency using NIS (Normalized Innovation Squared)

## Build and Run

```bash
mkdir build && cd build
cmake .. && make
./ukf_highway
```

## Project Structure

```text
src/
├── ukf.cpp & ukf.h          # UKF implementation
├── main.cpp                 # Highway simulation
├── highway.h                # Environment setup
└── tools.cpp & tools.h     # RMSE utilities

build/
├── final_rmse.txt           # Performance results
└── tuning_results.csv       # Parameter optimization log
```
