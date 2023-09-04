// See Chapter 2 of A Programmer's Introduction to Mathematics (https://pimbook.org)

// TODO:
//  1. Add option for random generation of Polynomial coefficients
//  2. Add option for random generation of test points
//  3. Rewrite this as a library with proper interface and structuring
//  4. Write tests
//  5. Add checks (e.g., interpolation points are actually different)
//  6. Make NewtonPolynomial::ddiff more efficient

#[derive(Debug, PartialEq)]
struct Point {
    x : f32,
    y : f32
}

trait PolyInterpolate<'a> {
    fn interpolate(points: &'a [Point]) -> Self;
}

trait PolyGetPoints {

    fn get_y(&self, x: &f32) -> f32;

    fn get_points(&self, xs : &[f32]) -> Vec<Point> {
        xs.iter()
          .map(|x| Point{x: *x, y: self.get_y(x)})
          .collect()
    }
}

#[derive(Debug)]
struct Polynomial<'a>(&'a[f32]);

impl<'a> Polynomial<'a> {
    fn new(coeffs : &'a[f32]) -> Self {
        Polynomial(coeffs)
    }
}

impl PolyGetPoints for Polynomial<'_> {
    fn get_y(&self, x: &f32) -> f32 {
        self.0.iter()
            .enumerate()
            .fold(0.0, |acc, c| acc + x.powi(c.0 as i32) * (c.1))
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
}

impl PolyGetPoints for LagrangePolynomial<'_> {
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
}

impl<'a> PolyInterpolate<'a> for LagrangePolynomial<'a> {
    fn interpolate(points: &'a [Point]) -> Self {
        LagrangePolynomial { bterms: LagrangePolynomial::get_bweights(points) }
    }
}

// See https://en.wikipedia.org/wiki/Newton_polynomial

#[derive(Debug)]
struct NewtonPolynomial<'a> {
    points : &'a [Point],
    ddiffs : Vec<f32> // divided differences
}

impl NewtonPolynomial<'_> {
    fn get_ddiffs(points : &[Point]) -> Vec<f32> {
        points.iter()
            .enumerate()
            .map(|j| Self::ddiff(0, j.0 as i32, points))
            .collect()
    }

    // compute divided differences with naive recursion
    fn ddiff(i: i32, j: i32, points: &[Point]) -> f32 {
        match (i - j).abs() {
            0 => points[i as usize].y,
            1 => (points[j as usize].y - points[i as usize].y) / (points[j as usize].x - points[i as usize].x),
            _ => (Self::ddiff(i + 1, j, points) - Self::ddiff(i, j - 1, points)) / (points[j as usize].x - points[i as usize].x)
        }
    }
}

impl<'a> PolyInterpolate<'a> for NewtonPolynomial<'a> {
    fn interpolate(points: &'a [Point]) -> Self {
        NewtonPolynomial {points : points, ddiffs : NewtonPolynomial::get_ddiffs(points)}
    }
}

impl PolyGetPoints for NewtonPolynomial<'_> {

    fn get_y(&self, x: &f32) -> f32 {
        if let Some(point) = self.points
            .iter()
            .find(|p| p.x == *x) {
                point.y
        } else {
            self.ddiffs[1..].iter()
                .zip(self.points.iter()
                        .map(|p| *x - p.x )
                        .scan(1_f32, |acc, x| {
                            *acc = *acc * x;
                            Some(*acc)
                        }) // Newton basis polynomials
                )
                .map(|z| *z.0 * z.1)
                .sum::<f32>()
            + self.ddiffs[0] // first divided difference in the sum doesn't have a multiplier
        }
    }
}


fn main() {
    let p: Polynomial = Polynomial::new(&[1.9, 9.2, 7.0]);
    let points = p.get_points(&[1.8,37.2,80.9]);

    let lp = LagrangePolynomial::interpolate(&points);
    let np = NewtonPolynomial::interpolate(&points);

    let test_points: Vec<f32> = (10u8..100u8).map(f32::from).collect();

    let count_lp = p.get_points(&test_points)
        .into_iter()
        .zip(lp.get_points(&test_points).into_iter())
        .filter(|x| (x.0.y - x.1.y).abs() > 0.01)
        .count();
    println!("{}", count_lp);

    let count_np = p.get_points(&test_points)
        .into_iter()
        .zip(np.get_points(&test_points).into_iter())
        .filter(|x| (x.0.y - x.1.y).abs() > 0.05)
        .count();
    println!("{}", count_np);
}