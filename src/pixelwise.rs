use std::f32::consts::PI;

use indicatif::{ProgressBar, ProgressStyle};
use mt19937::MT19937;
use nalgebra::DMatrix;
use rand_distr::{Distribution, Normal};

use crate::{
    rand::{
        cellseed, my_rand_poisson, myrand_gaussian_0_1, myrand_uniform_0_1, mysrand, NoisePrng,
    },
    utils::sq_dist,
};
use rayon::prelude::*;

fn render_pixel(
    img_in: &DMatrix<f32>,
    y_out: i32,
    x_out: i32,
    m_in: u32,
    n_in: u32,
    m_out: u32,
    n_out: u32,
    offset: u32,
    n_montecarlo: u32,
    grain_radius: f32,
    sigma_r: f32,
    sigma_filter: f32,
    x_a: f32,
    y_a: f32,
    x_b: f32,
    y_b: f32,
    lambdas: &[f32],
    exp_lambdas: &[f32],
    x_gaussians: &[f32],
    y_gaussians: &[f32],
) -> f32 {
    let normal_quantile = 3.0902; // Standard normal quantile for alpha=0.999
    let grain_radius_sq = grain_radius * grain_radius;
    let mut max_radius = grain_radius;
    let mut mu = 0.0;
    let mut sigma = 0.0;
    let sigma_sq;
    let mut curr_radius;
    let mut curr_radius_sq;

    let ag = 1.0 / (1.0 / grain_radius).ceil();
    let s_x = (n_out - 1) as f32 / (x_b - x_a);
    let s_y = (m_out - 1) as f32 / (y_b - y_a);

    // Pseudo-random number generators
    let mut p = NoisePrng::new();

    let mut pix_out = 0.0;

    // Convert from output grid (x_out, y_out) to input grid (x_in, y_in)
    let x_in = x_a + (x_out as f32 + 0.5) * (x_b - x_a) / n_out as f32;
    let y_in = y_a + (y_out as f32 + 0.5) * (y_b - y_a) / m_out as f32;

    // Calculate mu and sigma for the log-normal distribution
    if sigma_r > 0.0 {
        sigma = ((sigma_r / grain_radius).powi(2) + 1.0).ln().sqrt();
        sigma_sq = sigma.powi(2);
        mu = grain_radius.ln() - sigma_sq / 2.0;
        let log_normal_quantile = (mu + sigma * normal_quantile).exp();
        max_radius = log_normal_quantile;
    }

    // Monte Carlo simulation
    for i in 0..n_montecarlo {
        let x_gaussian = x_in + sigma_filter * x_gaussians[i as usize] / s_x;
        let y_gaussian = y_in + sigma_filter * y_gaussians[i as usize] / s_y;

        // Determine bounding box around the current shifted pixel
        let min_x = ((x_gaussian - max_radius) / ag).floor() as i32;
        let max_x = ((x_gaussian + max_radius) / ag).floor() as i32;
        let min_y = ((y_gaussian - max_radius) / ag).floor() as i32;
        let max_y = ((y_gaussian + max_radius) / ag).floor() as i32;

        // Iterate through cells

        'outer: for ncx in min_x..=max_x {
            for ncy in min_y..=max_y {
                let cell_corner_x = ag * ncx as f32;
                let cell_corner_y = ag * ncy as f32;

                // Seed the cell
                let seed = cellseed(ncx as u32, ncy as u32, offset);
                mysrand(&mut p, seed);

                // Compute Poisson parameters for the pixel that contains (x,y)
                let y_index = (cell_corner_y.max(0.0).min((m_in - 1) as f32)) as usize;
                let x_index = (cell_corner_x.max(0.0).min((n_in - 1) as f32)) as usize;
                let u = img_in[(y_index, x_index)];
                let u_ind = (u * (lambdas.len() as f32 - 1.0)).floor() as usize;
                let curr_lambda = lambdas[u_ind];
                let curr_exp_lambda = exp_lambdas[u_ind];

                // Draw number of points in the cell
                let n_cell = my_rand_poisson(&mut p, curr_lambda, curr_exp_lambda);

                // Process each point in the cell
                for _ in 0..n_cell {
                    let x_centre_grain = cell_corner_x + ag * myrand_uniform_0_1(&mut p);
                    let y_centre_grain = cell_corner_y + ag * myrand_uniform_0_1(&mut p);

                    if sigma_r > 0.0 {
                        // curr_radius =  fmin(exp(mu + sigma*myrand_gaussian_0_1(&p)),maxRadius); - c++
                        curr_radius = (mu + sigma * myrand_gaussian_0_1(&mut p))
                            .exp()
                            .min(max_radius);
                        curr_radius_sq = curr_radius.powi(2);
                    } else if sigma_r == 0.0 {
                        curr_radius_sq = grain_radius_sq;
                    } else {
                        panic!("sigma_r must be >= 0.0");
                    }

                    // test distance

                    if sq_dist(x_centre_grain, y_centre_grain, x_gaussian, y_gaussian)
                        < curr_radius_sq
                    {
                        pix_out += 1.0;
                        break 'outer;
                    }
                }
            }
        }
    }

    (pix_out / n_montecarlo as f32).clamp(0.0, 1.0)
}

#[derive(Debug, Clone)]
pub struct FilmGrainOptions {
    pub mu_r: f32,
    pub sigma_r: f32,
    pub sigma_filter: f32,
    pub n_montecarlo: usize,
    pub grain_seed: u32,
    pub m_out: usize,
    pub n_out: usize,
    pub x_a: f32,
    pub y_a: f32,
    pub x_b: f32,
    pub y_b: f32,
}

pub fn film_grain_rendering_pixel_wise(
    img_in: &DMatrix<f32>,
    film_grain_options: FilmGrainOptions,
) -> DMatrix<f32> {
    let grain_radius = film_grain_options.mu_r;
    let grain_std = film_grain_options.sigma_r;
    let sigma_filter = film_grain_options.sigma_filter;

    let n_montecarlo = film_grain_options.n_montecarlo;
    let mut rng = MT19937::default();

    // Draw random (Gaussian) translation vectors
    let normal_dist = Normal::new(0.0, sigma_filter as f64);

    let normal_dist = match normal_dist {
        Ok(n) => n,
        Err(e) => {
            println!("{}", e);
            panic!("erromax_radiusr!");
        }
    };

    let mut x_gaussian_list: Vec<f32> = vec![0.0; n_montecarlo];
    let mut y_gaussian_list: Vec<f32> = vec![0.0; n_montecarlo];

    for i in 0..n_montecarlo {
        x_gaussian_list[i] = normal_dist.sample(&mut rng) as f32;
        y_gaussian_list[i] = normal_dist.sample(&mut rng) as f32;
    }

    // Pre-calculate lambda and exp(-lambda) for each possible grey level
    const MAX_GREY_LEVEL: usize = 255;
    let epsilon_grey_level = 0.1;

    let mut lambda_list = vec![0.0; MAX_GREY_LEVEL + 1];
    let mut exp_lambda_list = vec![0.0; MAX_GREY_LEVEL + 1];

    let ag = 1.0 / (1.0 / grain_radius).ceil();
    for i in 0..=MAX_GREY_LEVEL {
        let u = (i as f32) / (MAX_GREY_LEVEL as f32 + epsilon_grey_level);
        let lambda_temp = -(ag * ag / (PI * (grain_radius * grain_radius + grain_std * grain_std)))
            * (1.0 - u).ln();
        lambda_list[i] = lambda_temp;
        exp_lambda_list[i] = (-lambda_temp).exp();
    }

    // Create output image
    let mut img_out =
        DMatrix::from_element(film_grain_options.m_out, film_grain_options.n_out, 0.0);

    // Get dimensions of img_out
    let nrows = img_out.nrows() as u32;
    let ncols = img_out.ncols() as u32;

    // Convert img_out to a mutable slice for parallel processing
    let img_out_slice = img_out.as_mut_slice();

    let total_pixels = nrows * ncols;
    let pb = ProgressBar::new(total_pixels as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{bar:70} {pos}/{len} ({percent}%)")
            .expect("Error setting progress bar style"),
    );

    img_out_slice
        .par_iter_mut()
        .enumerate()
        .for_each(|(index, pix_out)| {
            let i = index / ncols as usize; // Row index
            let j = index % ncols as usize; // Column index

            // Call render_pixel to compute the pixel value
            *pix_out = render_pixel(
                img_in,
                i as i32,
                j as i32,
                img_in.nrows() as u32,
                img_in.ncols() as u32,
                nrows,
                ncols,
                film_grain_options.grain_seed,
                n_montecarlo as u32,
                grain_radius,
                grain_std,
                sigma_filter,
                film_grain_options.x_a,
                film_grain_options.y_a,
                film_grain_options.x_b,
                film_grain_options.y_b,
                &lambda_list,
                &exp_lambda_list,
                &x_gaussian_list,
                &y_gaussian_list,
            );

            pb.inc(1);
        });

    // Return the modified img_out matrix
    img_out
}
