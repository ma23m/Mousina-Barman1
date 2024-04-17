// Define the differential equation dy/dt = f(t, y)
fn f(t: f64, y: f64) -> f64 {
    1.0 - t.powi(2) + y
}

// Fourth-order Runge-Kutta method
fn rk4(f: fn(f64, f64) -> f64, a: f64, b: f64, y0: f64, n: usize) -> (Vec<f64>, Vec<f64>) {
    let h = (b - a) / n as f64;
    let mut t = Vec::with_capacity(n + 1);
    let mut y = Vec::with_capacity(n + 1);
    t.push(a);
    y.push(y0);

    for i in 0..n {
        let k1 = h * f(t[i], y[i]);
        let k2 = h * f(t[i] + 0.5 * h, y[i] + 0.5 * k1);
        let k3 = h * f(t[i] + 0.5 * h, y[i] + 0.5 * k2);
        let k4 = h * f(t[i] + h, y[i] + k3);
        y.push(y[i] + (1.0 / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4));
        t.push(a + (i as f64 + 1.0) * h);
    }

    (t, y)
}

// Define the initial condition and interval
const A: f64 = 0.0;
const B: f64 = 2.0;
const Y0: f64 = 0.5;
const N: usize = 10;

fn main() {
    // Solve the differential equation using RK4
    let (t, y) = rk4(f, A, B, Y0, N);

    // Print the solution
    println!("t\t y");
    for i in 0..=N {
        println!("{:.2}\t {:.6}", t[i], y[i]);
    }
}
