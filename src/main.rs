mod pixelwise;
mod rand;
mod utils;

use crate::pixelwise::{film_grain_rendering_pixel_wise, FilmGrainOptions};
use image::{GenericImageView, ImageBuffer, Rgba};
use nalgebra::DMatrix;
use std::error::Error;

const MAX_COLOR_VALUE: u8 = 255;
const EPSILON_COLOR_LEVEL: f32 = 0.1;

fn main() -> Result<(), Box<dyn Error>> {
    // Load an image from a file
    let img_path = std::env::args().nth(1).ok_or("No file path provided")?;
    let output_path = std::env::args()
        .nth(2)
        .ok_or("No output file path provided")?;
    let img = image::open(&img_path)?;

    let normalization_factor = MAX_COLOR_VALUE as f32 + EPSILON_COLOR_LEVEL;

    // Convert the image to separate RGB matrices
    let (width, height) = img.dimensions();
    let mut img_matrix_r = DMatrix::from_element(height as usize, width as usize, 0.0);
    let mut img_matrix_g = DMatrix::from_element(height as usize, width as usize, 0.0);
    let mut img_matrix_b = DMatrix::from_element(height as usize, width as usize, 0.0);

    for (x, y, pixel) in img.pixels() {
        let red_value = pixel.0[0] as f32 / normalization_factor;
        let green_value = pixel.0[1] as f32 / normalization_factor;
        let blue_value = pixel.0[2] as f32 / normalization_factor;

        img_matrix_r[(y as usize, x as usize)] = red_value;
        img_matrix_g[(y as usize, x as usize)] = green_value;
        img_matrix_b[(y as usize, x as usize)] = blue_value;
    }

    // Print some statistics about the input matrices
    println!(
        "Input matrix (R) min: {}, max: {}",
        img_matrix_r.min(),
        img_matrix_r.max()
    );
    println!(
        "Input matrix (G) min: {}, max: {}",
        img_matrix_g.min(),
        img_matrix_g.max()
    );
    println!(
        "Input matrix (B) min: {}, max: {}",
        img_matrix_b.min(),
        img_matrix_b.max()
    );

    let x_a = 0.0;
    let y_a = 0.0;
    let x_b = width as f32;
    let y_b = height as f32;

    // Define film grain options
    let film_grain_options = FilmGrainOptions {
        mu_r: 0.1,
        sigma_r: 0.0,
        sigma_filter: 0.8,
        n_montecarlo: 800,
        grain_seed: 42,
        m_out: height as usize,
        n_out: width as usize,
        x_a,
        y_a,
        x_b,
        y_b,
    };

    // Perform film grain rendering for each channel
    let img_out_r = film_grain_rendering_pixel_wise(&img_matrix_r, film_grain_options.clone());
    let img_out_g = film_grain_rendering_pixel_wise(&img_matrix_g, film_grain_options.clone());
    let img_out_b = film_grain_rendering_pixel_wise(&img_matrix_b, film_grain_options);

    // Print some statistics about the output matrices
    println!(
        "Output matrix (R) min: {}, max: {}",
        img_out_r.min(),
        img_out_r.max()
    );
    println!(
        "Output matrix (G) min: {}, max: {}",
        img_out_g.min(),
        img_out_g.max()
    );
    println!(
        "Output matrix (B) min: {}, max: {}",
        img_out_b.min(),
        img_out_b.max()
    );

    // Create an output image buffer
    let mut output_image = ImageBuffer::new(width, height);

    // Convert the output DMatrix back to image
    for (index, ((&r_value, &g_value), &b_value)) in img_out_r
        .iter()
        .zip(img_out_g.iter())
        .zip(img_out_b.iter())
        .enumerate()
    {
        let row = index / img_out_r.ncols();
        let col = index % img_out_r.ncols();

        // Ensure the value is in the range [0, 1]
        let normalized_r = (r_value * normalization_factor).clamp(0.0, 255.0) as u8;
        let normalized_g = (g_value * normalization_factor).clamp(0.0, 255.0) as u8;
        let normalized_b = (b_value * normalization_factor).clamp(0.0, 255.0) as u8;

        let pixel = Rgba([normalized_r, normalized_g, normalized_b, 255]);

        output_image.put_pixel(col as u32, row as u32, pixel);
    }

    output_image.save(&output_path)?;

    println!("Image saved to {}", output_path);

    Ok(())
}
