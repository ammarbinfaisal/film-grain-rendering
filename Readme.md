# Realistic Film Grain Rendering

## Installation

### Prerequisites

- Rust (latest stable version)
- Cargo (comes with Rust)

### Steps

1. Clone this repository:
   ```
   git clone https://github.com/ammarbinfaisal/film-grain-rendering.git
   cd film-grain-rendering
   ```

2. Build the project:
   ```
   cargo build --release
   ```

## Usage

Run the program with the following command:

```
cargo run --release -- <input_image_path> <output_image_path>
```

Replace `<input_image_path>` with the path to your input image and `<output_image_path>` with the desired path for the output image.

## Configuration

You can adjust the film grain parameters in the `main.rs` file:

```rust
let film_grain_options = FilmGrainOptions {
    mu_r: 0.1,
    sigma_r: 0.0,
    sigma_filter: 0.8,
    n_montecarlo: 800,
    grain_seed: 42,
    // ... other options
};
```

- `mu_r`: Mean grain radius
- `sigma_r`: Standard deviation of grain radius
- `sigma_filter`: Standard deviation of the Gaussian filter
- `n_montecarlo`: Number of Monte Carlo iterations
- `grain_seed`: Seed for random number generation

## Project Structure

- `main.rs`: Entry point of the application
- `pixelwise.rs`: Contains the core film grain rendering algorithm
- `rand.rs`: Custom random number generation functions
- `utils.rs`: Utility functions

## Performance

The program uses parallel processing to improve performance. Progress is displayed using a progress bar.

## Dependencies

- `image`: For image processing
- `nalgebra`: For matrix operations
- `indicatif`: For progress bar
- `rayon`: For parallel processing

## Future Improvements

- Add command-line options for grain parameters
- Implement more grain patterns
- Add gainwise grain rendering
