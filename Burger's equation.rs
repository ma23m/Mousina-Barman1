use std::f64;

// Define parameters
const NU: f64 = 0.1;  // kinematic viscosity
const LENGTH: f64 = 2.0 * f64::consts::PI;  // spatial domain length
const NUM_PT: usize = 100;  // number of spatial points
const DX: f64 = LENGTH / NUM_PT as f64;  // spatial step size
const DT: f64 = 0.01;  // time step size

// Define Burger's equation function
fn burgers_equation(u: &Vec<f64>) -> Vec<f64> {
    let mut dudt = vec![0.0; NUM_PT];
    let du_dx = u.iter().skip(1).zip(u.iter().take(NUM_PT - 1)).map(|(x2, x1)| (x2 - x1) / DX).collect::<Vec<f64>>();
    for i in 1..NUM_PT - 1 {
        dudt[i] = -u[i] * du_dx[i - 1] - u[i] * (u[i + 1] - u[i - 1]) / (2.0 * DX) + NU * (u[i + 1] - 2.0 * u[i] + u[i - 1]) / DX.powi(2);
    }
    dudt
}

// Define initial condition (cosine function)
fn initial_condition(x: f64) -> f64 {
    f64::cos(x)
}

// Fourth-order Runge-Kutta method
fn rk4_step(u: &Vec<f64>) -> Vec<f64> {
    let k1 = burgers_equation(u);
    let k2 = burgers_equation(&u.iter().zip(&k1).map(|(u_val, k_val)| u_val + 0.5 * DT * k_val).collect());
    let k3 = burgers_equation(&u.iter().zip(&k2).map(|(u_val, k_val)| u_val + 0.5 * DT * k_val).collect());
    let k4 = burgers_equation(&u.iter().zip(&k3).map(|(u_val, k_val)| u_val + DT * k_val).collect());
    u.iter().zip(&k1).zip(&k2).zip(&k3).zip(&k4)
        .map(|((((u_val, k1_val), k2_val), k3_val), k4_val)| u_val + (DT / 6.0) * (k1_val + 2.0 * k2_val + 2.0 * k3_val + k4_val))
        .collect()
}

// Main function
fn main() {
    // Initialize arrays
    let mut u = vec![0.0; NUM_PT];
    let mut x = vec![0.0; NUM_PT];
    for i in 0..NUM_PT {
        x[i] = i as f64 * DX;
        u[i] = initial_condition(x[i]);
    }

    // Print the header for the table
    println!("Step\t{}", (0..10).map(|i| format!("u[{}]", i)).collect::<Vec<String>>().join("\t"));

    // Main loop using RK4
    for step in 0..10 {
        // Print the solution at the current step
        print!("{}\t", step);
        for i in 0..10 {
            print!("{:.6}\t", u[i]);
        }
        println!();
        
        // Perform RK4 step
        u = rk4_step(&u);
    }
}
