// See Chapter 2 of A Programmer's Introduction to Mathematics (https://pimbook.org)

// TODO:
//  1. Add option for random generation of Polynomial coefficients
//  2. Add option for random generation of test points
//  3. Rewrite this as a library with proper interface and structuring
//  4. Write tests
//  5. Add checks (e.g., interpolation points are actually different)
//  6. Handle case where interpolating point is fed into get_points()

use std::vec;

#[derive(Debug, PartialEq)]
struct Point {
    x : f32,
    y : f32
}

trait PolyInterpolate {
    fn interpolate(points: &[Point]) -> Self; 
}

#[derive(Debug)]
struct Polynomial(Vec<f32>);

impl Polynomial {
    fn new(coeffs : Vec<f32>) -> Self {
        Polynomial(coeffs)
    }

    fn get_y(&self, x: &f32) -> f32 {
        self.0.iter()
            .enumerate()
            .fold(0.0, |acc, c| acc + x.powi(c.0 as i32) * (c.1))
    }

    fn get_points(&self, xs : &[f32]) -> Vec<Point> {
        xs.iter()
          .map(|x| Point{x: *x, y: self.get_y(x)})
          .collect()
    }
}

// True Barycentric form is used for computing the Lagrange interpolating polynomial.
// See https://en.wikipedia.org/wiki/Lagrange_polynomial#Barycentric_form

#[derive(Debug)]
struct Bterm {
    w : f32,
    p : Point
}

#[derive(Debug)]
struct LagrangePolynomial {
    bterms : Vec<Bterm>
}

impl LagrangePolynomial {
    fn get_bweights(points: &[Point]) -> Vec<Bterm> {
        points.iter()
                .map(|point: &Point|
                    Bterm {
                        p: Point {x: point.x, y: point.y},
                        w: points
                            .iter()
                            .map(|p: &Point | point.x - p.x )
                            .filter(|p| *p != 0_f32)
                            .product::<f32>()
                    }
                )
                .collect()
    }

    fn get_y(&self, x:&f32) -> f32 {
        let terms: (f32, f32) = self.bterms
            .iter()
            .fold((0.0, 0.0),
                |acc, bterm| {
                    let temp = (x - bterm.p.x) * bterm.w;
                    (acc.0 + (bterm.p.y / temp), acc.1 + (1.0 / temp))
                }
            );

        terms.0 / terms.1
    }

    fn get_points(&self, xs : &[f32]) -> Vec<Point> {
        xs.iter()
          .map(|x| Point {x: *x, y: self.get_y(x)})
          .collect()
    }
}

impl PolyInterpolate for LagrangePolynomial {
    fn interpolate(points: &[Point]) -> Self {
        LagrangePolynomial { bterms: LagrangePolynomial::get_bweights(points) }
    }
}


// See https://en.wikipedia.org/wiki/Newton_polynomial

#[derive(Debug)]
struct NewtonPolynomial {}

impl PolyInterpolate for NewtonPolynomial {
    fn interpolate(_points: &[Point]) -> Self {
        NewtonPolynomial {}
    }
}

fn main() {
    let p: Polynomial = Polynomial::new(vec!(1.9, 9.2, 7.0));
    let points = p.get_points(&[1.8,37.2,80.9]);

    let lp = LagrangePolynomial::interpolate(&points);

    let test_points: Vec<f32> = (10u8..100u8).map(f32::from).collect();
    let count = p.get_points(&test_points)
        .into_iter()
        .zip(lp.get_points(&test_points).into_iter())
        .filter(|x| (x.0.y - x.1.y).abs() > 0.01)
        .count();
    println!("{}", count);

    println!("{:?}", lp.get_points(&[44_f32]));
}