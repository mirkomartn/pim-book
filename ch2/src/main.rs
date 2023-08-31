// See Chapter 2 of A Programmer's Introduction to Mathematics (https://pimbook.org)

// TODO:
//  1. Add option for random generation of Polynomial coefficients
//  2. Add option for random generation of test points
//  3. Rewrite this as a library with proper interface and structuring
//  4. Write tests
//  5. Add checks (e.g., interpolation points are actually different)

#[derive(Debug, PartialEq)]
struct Point {
    x : f32,
    y : f32
}

trait PolyInterpolate<'a> {
    fn interpolate(points: &'a [Point]) -> Self;
}

#[derive(Debug)]
struct Polynomial<'a>(&'a[f32]);

impl<'a> Polynomial<'a> {
    fn new(coeffs : &'a[f32]) -> Self {
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
struct Bterm<'a> {
    w : f32,
    p : &'a Point
}

#[derive(Debug)]
struct LagrangePolynomial<'a> {
    bterms : Vec<Bterm<'a>>
}

impl<'a> LagrangePolynomial<'a> {
    fn get_bweights(points: &'a [Point]) -> Vec<Bterm> {
        points.iter()
                .map(|point: &Point|
                    Bterm {
                        p: point,
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
        // check if this is one of the interpolation points
        if let Some(bweight) = self.bterms
            .iter()
            .find(|b| b.p.x == *x) {
                bweight.p.y
        }
        // else compute y
        else {
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
    }

    fn get_points(&self, xs : &[f32]) -> Vec<Point> {
        xs.iter()
          .map(|x| Point {x: *x, y: self.get_y(x)})
          .collect()
    }
}

impl<'a> PolyInterpolate<'a> for LagrangePolynomial<'a> {
    fn interpolate(points: &'a [Point]) -> Self {
        LagrangePolynomial { bterms: LagrangePolynomial::get_bweights(points) }
    }
}


// See https://en.wikipedia.org/wiki/Newton_polynomial

#[derive(Debug)]
struct NewtonPolynomial {}

// impl<'a> PolyInterpolate<'a> for NewtonPolynomial<'a> {
//     fn interpolate(_points: &'a [Point]) -> Self {
//         NewtonPolynomial {}
//     }
// }


fn main() {
    let p: Polynomial = Polynomial::new(&[1.9, 9.2, 7.0]);
    let points = p.get_points(&[1.8,37.2,80.9]);

    let lp = LagrangePolynomial::interpolate(&points);

    let test_points: Vec<f32> = (10u8..100u8).map(f32::from).collect();
    let count = p.get_points(&test_points)
        .into_iter()
        .zip(lp.get_points(&test_points).into_iter())
        .filter(|x| (x.0.y - x.1.y).abs() > 0.01)
        .count();
    println!("{}", count);
}