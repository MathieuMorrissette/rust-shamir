extern crate num;
extern crate num_bigint;

use num_bigint::BigInt;
use num::integer::Integer;
use num_bigint::ToBigInt;


fn _eval_at(poly: Vec<BigInt>, x: i32, prime: BigInt) -> BigInt{
    let mut accum : BigInt = 0.into();

    for coeff in poly.iter().rev() {
        accum = accum * x;
        accum = accum + coeff;
        accum = accum % &prime;
    }

    accum
}

fn _extended_gcd(mut a: BigInt,mut b: BigInt) -> (BigInt, BigInt) {
    let mut x:BigInt = 0.into();
    let mut last_x:BigInt = 1.into();
    let mut y:BigInt = 1.into();
    let mut last_y:BigInt = 0.into();

    let zero:BigInt = 0.into();

    while b != zero {
        let quot = a.div_floor(&b);
        
        let temp1 = a.mod_floor(&b);
        let temp2 = b.clone();
        
        a = temp2;
        b = temp1;        

        let temp3 = last_x - &quot * &x;        
        let temp4 = x.clone();

        x = temp3;
        last_x = temp4;

        let temp5 = last_y - &quot * &y;        
        let temp6 = y.clone();

        y = temp5;
        last_y = temp6;
    }

    (last_x, last_y)
}

fn _divmod(num: BigInt, den:BigInt, p:BigInt) -> BigInt {

    let (inv, _) = _extended_gcd(den,p);

    num * inv
}


fn _lagrange_interpolate(x: BigInt, x_s: Vec<BigInt> , y_s: Vec<BigInt>, p: BigInt) -> BigInt {
    let k = x_s.len();

    let mut nums: Vec<BigInt> = Vec::new();
    let mut dens: Vec<BigInt> = Vec::new();    

    for i in 0..k{
       // missing stuff here
        let mut others = x_s.clone();

        let cur = others.remove(i);

        println!("{}", cur);
        
        for o in others{
            nums.push(pi(vec![&x - &o]));
            dens.push(pi(vec![&cur - &o]));
        }        
    }

    let mut den: BigInt = pi(dens.clone());

    println!("{}", den);

    let mut num : BigInt = 0.into();

    for i in 0..k{
        num = num + _divmod(&nums[i] * &den * &y_s[i] % &p, dens[i].clone(), p.clone());
    }

    println!("{}", num);

    (_divmod(num.clone(), den.clone(), p.clone()) + &p) % &p
}

fn pi(vals: Vec<BigInt>) -> BigInt{
    let mut accum: BigInt = 1.into();
    for v in vals {
        accum = accum * v;
    }

    accum
}

// 12th Mersenne Prime
const _PRIME:i128 = 170141183460469231731687303715884105727;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_eval() {
        
        let mut poly: Vec<BigInt> = Vec::new();        
        
        poly.push(146954493737107699617854775963500452322_i128.into());
        poly.push(19758689580932147179251136439148532855_i128.into());
        poly.push(49127370899614395569901674483194295184_i128.into());

        
        let res = _eval_at(poly, 2,  170141183460469231731687303715884105727_i128.into());
        
        println!("{}", res);
        assert_eq!(res, 42698989576491112792589139342806487314_i128.into());         
    }

    #[test]
    fn test_gcd(){
        let tuple = _extended_gcd(2.into(),170141183460469231731687303715884105727_i128.into());

        println!("{}", tuple.0);
        println!("{}", tuple.1);

        assert_eq!(tuple.0, -85070591730234615865843651857942052863_i128.to_bigint().unwrap());
        assert_eq!(tuple.1, 1.into());

    }

    #[test]
    fn test_divmod()
    {
        let num:BigInt = 137288411751325586557619846070933607306_i128.into();
        let den:BigInt = 2.into();
        
        let res = _divmod(num, den, 170141183460469231731687303715884105727_i128.into());

        assert_eq!(res, BigInt::parse_bytes(b"-11679206425389363299939391086680605274249979452515224891674157973087535017078", 10).unwrap());
    }

    #[test]
    fn test_pi(){
        let mut vals: Vec<BigInt> = Vec::new();

        vals.push(BigInt::parse_bytes(b"-2", 10).unwrap());
        vals.push(BigInt::parse_bytes(b"-3", 10).unwrap());

        let res = pi(vals);

        assert_eq!(res, 6.into());

    }

    #[test]
    fn test_lagrange(){
        let x: BigInt = 0.into();

        let mut x_s: Vec<BigInt> = Vec::new();

        x_s.push(1.into());
        x_s.push(2.into());
        x_s.push(3.into());

        let mut y_s: Vec<BigInt> = Vec::new();

        y_s.push(22873837885008894210735194230010426715_i128.into());
        y_s.push(20112177421602409233155423655360946858_i128.into());
        y_s.push(102484273664769396574413799315550648172_i128.into());

        let res = _lagrange_interpolate(x, x_s, y_s, 170141183460469231731687303715884105727_i128.into());

        println!("{}", res);
        assert_eq!(res, 110769255054988851507153111039499087743_i128.into());

    }

}